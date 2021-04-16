# Memory management in acados:
Guidelines on how memory should be assigned for an `acados` structure `astruct`.

There are two functions: `astruct_calculate_size()`, `astruct_assign()`.

## `astruct_calculate_size()`
Must return a multiple of 8 to keep the pointer aligned to 8 bytes when allocating substructures.
Thus, it should end with:
```
    make_int_multiple_of(8, &size);
```


## `astruct_assign()`
Should assign its members in the following order:

- Align to 8 bytes, i.e.:
```
    align_char_to(8, &c_ptr);
```

- Assign structure itself, i.e.:
```
    astruct *as_instance = (astruct *) c_ptr;
    c_ptr += sizeof(astruct);
```

- Assign pointers to substructures, if `astruct` contains an array of `bstruct`s, e.g.
```
    mem->bstructs = (void **) c_ptr;
    c_ptr += N*sizeof(void *);
```

- Align to 8 bytes, since `astruct` might contain `int`s and the pointers were assigned.


- Assign "substructures", i.e. structures that `astruct` has pointers to:
```
    // assume astruct contains an instance of bstruct
    as_instance->bstruct = bstruct_assign(c_ptr,...);
    c_ptr += bstruct_calculate_size(astruct);
```

- Note:
    - since calculate_size returns multiple of 8, c_ptr is still aligned to 8 bytes.
    - `blasfeo_dmat_structs`, `blasfeo_dvec_structs` can be seen as "substructures" here.


- Assign doubles (are 8 bytes anyway)
```
    assign_and_advance_double(n_doubles, &as_instance->doubles, &c_ptr);
```

- Assign pointers (4 bytes on 32 Bit)
```
    assign_and_advance_double_ptrs(n_pointers, &as_instance->double_pointer, &c_ptr);
```

- Align to 8 bytes, can be skipped if no pointers have been assigned


- Assign integers
```
    assign_and_advance_int(n_integers, &as_instance->ints, &c_ptr);
```

- Align to 64 bytes, i.e.:
```
    align_char_to(64, &c_ptr);
```

- Assign `blasfeo_dmat_mem` (are multiple of 64 Bytes)
```
    assign_and_advance_blasfeo_dmat_mem(nrow, ncol, &as_instance->blasfeo_mat, &c_ptr);
```

- Assign `blasfeo_vec_mem` (are multiple of 32 Bytes)
```
    assign_and_advance_blasfeo_dvec_mem(n, &as_instance->blasfeo_vec, &c_ptr);
```

- Align c_ptr to 8 byte here for nested assigns, see "substructures"
   - relevant if no blasfeo_mems are in `astruct`