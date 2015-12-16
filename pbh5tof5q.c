/*
 * 
 * Copyright (c) 2011, Jue Ruan <ruanjue@gmail.com>
 *
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <hdf5.h>
#include "list.h"
#include "file_reader.h"
#include "hashset.h"

int usage(){
	printf(
	"Convert pacbio hdf5 reads file into fast5q format\n"
	"'fast5q' inherits most features of fastq, the only difference\n"
	"   is that the quality line contain five times quality values\n"
	"   Suppose read length is L, its quality line likes:\n"
	"   L QualityValues fellowing by L SubstitutionQV, L InsertionQV, L DeletionQV, L MergeQV, L SubstitutionTag, and L DeletionTag\n"
	"   Quality value is encoded by phred score + 33, region from '!'(0) to '{'(90)\n"
	"   The original head of deletionQV is set to '~'(93+33)\n"
	"Author: Jue Ruan <ruanjue@gmail.com>\n"
	"Version: 1.0, 2015-06-13\n"
	"Usage: pbh5tof5q [options] <pacbio_h5_bas_file> >output.f5q\n"
	"Options:\n"
	" -l <int>    Min length of subreads, [100]\n"
	" -q <float>  Min read quality (RQ), [0.6]\n"
	" -b <int>    Batch size of loadinging reads, [1000]\n"
	//" -L          Output the longest subread of a hole\n"
	" -h          Show this document\n"
	"\n"
	);
	return 1;
}

typedef struct {
	uint32_t hole;
	uint32_t length;
	uint32_t x, y;
	int rq;
} pb_read_t;
define_list(pbreadv, pb_read_t);

static const char* qvspace[] = {
	"/PulseData/BaseCalls/QualityValue",
	"/PulseData/BaseCalls/SubstitutionQV",
	"/PulseData/BaseCalls/InsertionQV",
	"/PulseData/BaseCalls/DeletionQV",
	"/PulseData/BaseCalls/MergeQV",
	"/PulseData/BaseCalls/SubstitutionTag",
	"/PulseData/BaseCalls/DeletionTag"
};

void pbh5_pulse_region_to_pbregs(b32list *h5regs, uuhash *hash, pbreadv *rds){
	u32list *hqs;
	uint64_t i;
	uint32_t hole, idx, hx, hy, j, x, y, b, e;
	int rq;
	hole = hx = hy = x = y = rq = 0;
	hqs = init_u32list(8);
	for(i=0;i<=h5regs->size;i+=5){
		if(i == h5regs->size || h5regs->buffer[i] != (int)hole){
			if(hqs->size){
				if((idx = kv_get_uuhash(hash, hole)) != 0xFFFFFFFFU){
					for(j=0;j<hqs->size;j+=2){
						b = hqs->buffer[j+0] < hx? hx : hqs->buffer[j+0];
						e = hqs->buffer[j+1] > hy? hy : hqs->buffer[j+1];
						if(b >= e) continue;
						if(b + (y - x) < e){
							x = b; y = e;
						}
					}
					rds->buffer[idx].x  = x;
					rds->buffer[idx].y  = y;
					rds->buffer[idx].rq = rq;
				}
			}
			if(i == h5regs->size) break;
			hole = h5regs->buffer[i];
			hx = hy = x = y = rq = 0;
			clear_u32list(hqs);
		}
		if(h5regs->buffer[i+1] == 2){
			hx = h5regs->buffer[i+2];
			hy = h5regs->buffer[i+3];
			rq = h5regs->buffer[i+4];
		} else if(h5regs->buffer[i+1] == 1){
			push_u32list(hqs, h5regs->buffer[i+2]);
			push_u32list(hqs, h5regs->buffer[i+3]);
		}
	}
	free_u32list(hqs);
}

int main(int argc, char **argv){
	pbreadv *rds;
	pb_read_t *rd;
	uuhash *hole2idx;
	char *bas_file, *bax_file, *chp;
	String *libname, *filepath;
	cplist *parts;
	hid_t file, dataset, dataspace, datatype, mem_dt, mem_ds;
	hsize_t dims[1], hi, n_part, dsize;
	hsize_t cnt[1], off[2];
	u32list *pbidxs;
	u32list *pblens;
	b32list *h5regs;
	u8list  *bases, *qvs[7];
	uint64_t offset, length, boff;
	uint32_t i, j, idx, beg, end, batch_size, min_rd, longest;
	int c, file_idx, min_rq;
	batch_size = 1000;
	min_rd = 100;
	min_rq = 600;
	longest = 0;
	while((c = getopt(argc, argv, "hLl:q:b:")) != -1){
		switch(c){
			case 'h': return usage();
			case 'L': longest = 1; break;
			case 'l': min_rd = atoi(optarg); break;
			case 'q': min_rq = atof(optarg) * 1000; break;
			case 'b': batch_size = atoi(optarg); break;
			default: return usage();
		}
	}
	if(argc == optind) return usage();
	if(batch_size == 0) batch_size = 1;
	libname  = init_string(16);
	filepath = init_string(16);
	parts = init_cplist(16);
	for(file_idx=optind;file_idx<argc;file_idx++){
		bas_file = argv[file_idx];
		if(!file_exists(bas_file)){
			fprintf(stderr, " -- BAS File '%s' doesn't exist in %s -- %s:%d --\n", bas_file, __FUNCTION__, __FILE__, __LINE__);
			continue;
		}
		clear_string(libname);
		clear_string(filepath);
		if((chp = rindex(bas_file, '/'))){
			append_string(filepath, bas_file, chp + 1 - bas_file);
		}
		if(strcmp(".bas.h5", bas_file + strlen(bas_file) - 7) == 0){
			append_string(libname, bas_file + filepath->size, strlen(bas_file) - filepath->size - 7);
		} else {
			append_string(libname, bas_file + filepath->size, strlen(bas_file) - filepath->size);
		}
		file = H5Fopen((const char*)bas_file, H5F_ACC_RDONLY, H5P_DEFAULT);
		dataset = H5Dopen(file, "/MultiPart/Parts", H5P_DEFAULT);
		dataspace = H5Dget_space(dataset);
		H5Sget_simple_extent_dims(dataspace, dims, NULL);
		n_part = H5Sget_simple_extent_npoints(dataspace);
		clear_cplist(parts); encap_cplist(parts, n_part); parts->size = n_part;
		datatype = H5Dget_type(dataset);
		//H5Dvlen_get_buf_size(dataset, datatype, dataspace, &vsize);
		mem_dt = H5Tcopy(datatype);
		H5Dread(dataset, mem_dt, H5S_ALL, H5S_ALL, H5P_DEFAULT, parts->buffer);
		H5Tclose(mem_dt);
		H5Tclose(datatype);
		H5Sclose(dataspace);
		H5Dclose(dataset);
		H5Fclose(file);
		for(hi=0;hi<n_part;hi++){
			bax_file = catstr(2, filepath->string, parts->buffer[hi]);
			file = H5Fopen((const char*)bax_file, H5F_ACC_RDONLY, H5P_DEFAULT);
			if(file < 0){
				fprintf(stderr, " -- BAX File '%s' doesn't exist in %s -- %s:%d --\n", bax_file, __FUNCTION__, __FILE__, __LINE__);
				continue;
			}

			dataset = H5Dopen(file, "/PulseData/BaseCalls/ZMW/HoleNumber", H5P_DEFAULT);
			dataspace = H5Dget_space(dataset);
			dsize = H5Sget_simple_extent_npoints(dataspace);
			datatype = H5Dget_type(dataset);
			mem_dt = H5Tcopy(datatype);
			pbidxs = init_u32list(dsize);
			H5Dread(dataset, mem_dt, H5S_ALL, H5S_ALL, H5P_DEFAULT, pbidxs->buffer);
			pbidxs->size = dsize;
			H5Tclose(mem_dt);
			H5Tclose(datatype);
			H5Sclose(dataspace);
			H5Dclose(dataset);
			rds = init_pbreadv(pbidxs->size);
			hole2idx = init_uuhash(1023);
			for(i=0;i<pbidxs->size;i++){
				push_pbreadv(rds, (pb_read_t){pbidxs->buffer[i], 0, 0, 0, 0});
				kv_put_uuhash(hole2idx, pbidxs->buffer[i], i);
			}
			free_u32list(pbidxs);

			dataset = H5Dopen(file, "/PulseData/BaseCalls/ZMW/NumEvent", H5P_DEFAULT);
			dataspace = H5Dget_space(dataset);
			dsize = H5Sget_simple_extent_npoints(dataspace);
			datatype = H5Dget_type(dataset);
			mem_dt = H5Tcopy(datatype);
			pblens = init_u32list(dsize);
			H5Dread(dataset, mem_dt, H5S_ALL, H5S_ALL, H5P_DEFAULT, pblens->buffer);
			pblens->size = dsize;
			H5Tclose(mem_dt);
			H5Tclose(datatype);
			H5Sclose(dataspace);
			H5Dclose(dataset);
			for(i=0;i<pblens->size;i++) rds->buffer[i].length = pblens->buffer[i];
			free_u32list(pblens);

			dataset = H5Dopen(file, "/PulseData/Regions", H5P_DEFAULT);
			dataspace = H5Dget_space(dataset);
			dsize = H5Sget_simple_extent_npoints(dataspace);
			datatype = H5Dget_type(dataset);
			mem_dt = H5Tcopy(datatype);
			h5regs = init_b32list(dsize);
			H5Dread(dataset, mem_dt, H5S_ALL, H5S_ALL, H5P_DEFAULT, h5regs->buffer);
			h5regs->size = dsize;
			H5Tclose(mem_dt);
			H5Tclose(datatype);
			H5Sclose(dataspace);
			H5Dclose(dataset);
			pbh5_pulse_region_to_pbregs(h5regs, hole2idx, rds);
			free_b32list(h5regs);
			free_uuhash(hole2idx);

			offset = 0;
			for(idx=0;idx<rds->size;idx+=batch_size){
				beg = idx; end = idx + batch_size; if(end > rds->size) end = rds->size;
				length = 0;
				for(i=beg;i<end;i++) length += rds->buffer[i].length;
				bases  = init_u8list(length);
				qvs[0] = init_u8list(length);
				qvs[1] = init_u8list(length);
				qvs[2] = init_u8list(length);
				qvs[3] = init_u8list(length);
				qvs[4] = init_u8list(length);
				qvs[5] = init_u8list(length);
				qvs[6] = init_u8list(length);
				off[0] = offset;
				off[1] = 0;
				cnt[0] = length;
				{
					dataset = H5Dopen(file, "/PulseData/BaseCalls/Basecall", H5P_DEFAULT);
					dataspace = H5Dget_space(dataset);
					datatype = H5Dget_type(dataset);
					mem_dt = H5Tcopy(datatype);
					mem_ds = H5Scopy(dataspace);
					H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, off, NULL, cnt, NULL);
					H5Sselect_hyperslab(mem_ds, H5S_SELECT_SET, off + 1, NULL, cnt, NULL);
					H5Dread(dataset, mem_dt, mem_ds, dataspace, H5P_DEFAULT, bases->buffer);
					bases->size = length;
					H5Tclose(mem_dt);
					H5Tclose(datatype);
					H5Sclose(mem_ds);
					H5Sclose(dataspace);
					H5Dclose(dataset);
				}
				for(i=0;i<7;i++){
					dataset = H5Dopen(file, qvspace[i], H5P_DEFAULT);
					dataspace = H5Dget_space(dataset);
					datatype = H5Dget_type(dataset);
					mem_dt = H5Tcopy(datatype);
					mem_ds = H5Scopy(dataspace);
					H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, off, NULL, cnt, NULL);
					H5Sselect_hyperslab(mem_ds, H5S_SELECT_SET, off + 1, NULL, cnt, NULL);
					H5Dread(dataset, mem_dt, mem_ds, dataspace, H5P_DEFAULT, qvs[i]->buffer);
					qvs[i]->size = length;
					if(i < 5){
						for(j=0;j<length;j++){
							if(qvs[i]->buffer[j] >= 93) qvs[i]->buffer[j] = 93 + 33;
							else qvs[i]->buffer[j] += 33;
						}
					}
					H5Tclose(mem_dt);
					H5Tclose(datatype);
					H5Sclose(mem_ds);
					H5Sclose(dataspace);
					H5Dclose(dataset);
				}
				offset += length;
				boff = 0;
				for(i=beg;i<end;i++){
					rd = ref_pbreadv(rds, i);
					if(rd->length < min_rd || rd->rq < min_rq){
					} else {
						printf("@%s/%d/%d_%d RQ=%0.3f\n", libname->buffer, rd->hole, rd->x, rd->y, rd->rq / 1000.0);
						for(j=rd->x;j<rd->y;j++) putc(bases->buffer[boff + j], stdout);
						printf("\n+\n");
						for(j=rd->x;j<rd->y;j++) putc(qvs[0]->buffer[boff + j], stdout);
						for(j=rd->x;j<rd->y;j++) putc(qvs[1]->buffer[boff + j], stdout);
						for(j=rd->x;j<rd->y;j++) putc(qvs[2]->buffer[boff + j], stdout);
						qvs[3]->buffer[boff + 0] = 93 + 33; // '~'
						for(j=rd->x;j<rd->y;j++) putc(qvs[3]->buffer[boff + j], stdout);
						for(j=rd->x;j<rd->y;j++) putc(qvs[4]->buffer[boff + j], stdout);
						for(j=rd->x;j<rd->y;j++) putc(qvs[5]->buffer[boff + j], stdout);
						for(j=rd->x;j<rd->y;j++) putc(qvs[6]->buffer[boff + j], stdout);
						putc('\n', stdout);
					}
					boff += rd->length;
				}
				free_u8list(bases);
				free_u8list(qvs[0]);
				free_u8list(qvs[1]);
				free_u8list(qvs[2]);
				free_u8list(qvs[3]);
				free_u8list(qvs[4]);
				free_u8list(qvs[5]);
				free_u8list(qvs[6]);
			}
			free_pbreadv(rds);
			H5Fclose(file);
			free(parts->buffer[hi]);
			free(bax_file);
		}
	}
	free_string(libname);
	free_string(filepath);
	free_cplist(parts);
	return 0;
}
