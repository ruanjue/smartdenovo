#include "mem_share.h"

typedef struct {
	char *str;
	int val;
} Type1;

size_t type1_count(void *obj, int idx){ if(idx == 0){ return strlen(((Type1*)obj)->str) + 1; } else return 0; }
//const obj_desc_t type1_obj_desc = {sizeof(Type1), 1, {1}, {offsetof(Type1, str)}, {&OBJ_DESC_DATA}, type1_count};
const obj_desc_t type1_obj_desc = {sizeof(Type1), 1, {1}, {offsetof(Type1, str)}, {&OBJ_DESC_CHAR_ARRAY}, NULL};

typedef struct {
	int a, b, c;
	Type1 d1, d2[10], *d3, *d4[10], **d5;
	char **strs;
	int d3len, d5len;
} Type2;

size_t type2_count(void *obj, int idx){
	switch(idx){
		case 0: return 1;
		case 1: return 10;
		case 2: return ((Type2*)obj)->d3len;
		case 3: return 10;
		case 4: return ((Type2*)obj)->d5len;
		default: return 10;
	}
}

const obj_desc_t type2_obj_desc = {sizeof(Type2), 6, {0, 0, 1, 2, 3, 3}, {offsetof(Type2, d1), offsetof(Type2, d2), offsetof(Type2, d3), offsetof(Type2, d4), offsetof(Type2, d5), offsetof(Type2, strs)}, {&type1_obj_desc, &type1_obj_desc, &type1_obj_desc, &type1_obj_desc, &type1_obj_desc, &OBJ_DESC_CHAR_ARRAY}, type2_count};

int main(){
	Type2 *t2, *t3;
	t2 = calloc(1, sizeof(Type2));
	int idx = 0;
	t2->d1.val = idx ++;
	t2->d1.str = strdup("d1");
	int i;
	for(i=0;i<10;i++){
		t2->d2[i].val = idx ++;
		t2->d2[i].str = strdup("d2");
	}
	t2->d3len = 10;
	t2->d3 = malloc(sizeof(Type1) * t2->d3len);
	for(i=0;i<t2->d3len;i++){
		t2->d3[i].val = idx ++;
		t2->d3[i].str = strdup("d3");
	}
	for(i=0;i<10;i++){
		t2->d4[i] = malloc(sizeof(Type1));
		t2->d4[i]->val = idx ++;
		t2->d4[i]->str = strdup("d4");
	}
	t2->d5len = 10;
	t2->d5 = malloc(sizeof(Type1*) * t2->d5len);
	for(i=0;i<t2->d5len;i++){
		t2->d5[i] = malloc(sizeof(Type1));
		t2->d5[i]->val = idx ++;
		t2->d5[i]->str = strdup("d5");
	}
	t2->strs = malloc(sizeof(char*) * 10);
	for(i=0;i<10;i++){
		t2->strs[i] = malloc(32);
		sprintf(t2->strs[i], "strs[%d,%d]", i, idx ++);
	}
	size_t aux_data, size, cnt, mem_type;
	FILE *file;
	size = mem_size_obj(t2, 1, &type2_obj_desc, 0, 1);
	fprintf(stdout, " -- size = %d in %s -- %s:%d --\n", (int)size, __FUNCTION__, __FILE__, __LINE__);
	aux_data = 1000999900;
	file = fopen("test.mem_share", "w");
	size = mem_dump_free_obj_file(t2, 1, &type2_obj_desc, 1, aux_data, file);
	fclose(file);
	fprintf(stdout, " -- size = %d in %s -- %s:%d --\n", (int)(size - 4 * sizeof(size_t)), __FUNCTION__, __FILE__, __LINE__);
	t3 = mem_read_obj_file(&type2_obj_desc, "test.mem_share", &mem_type, &cnt, &aux_data);
	fprintf(stdout, " -- aux_data = %d in %s -- %s:%d --\n", (int)aux_data, __FUNCTION__, __FILE__, __LINE__);
	printf("%d %s\n", t3->d1.val, t3->d1.str);
	for(i=0;i<10;i++) printf("%d %s\n", t3->d2[i].val, t3->d2[i].str);
	for(i=0;i<t3->d3len;i++) printf("%d %s\n", t3->d3[i].val, t3->d3[i].str);
	for(i=0;i<10;i++) printf("%d %s\n", t3->d4[i]->val, t3->d4[i]->str);
	for(i=0;i<t3->d5len;i++) printf("%d %s\n", t3->d5[i]->val, t3->d5[i]->str);
	for(i=0;i<10;i++) printf("%s\n", t3->strs[i]);
	free(t3);
	return 0;
}
