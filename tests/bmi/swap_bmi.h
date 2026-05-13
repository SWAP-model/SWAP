/* SS-BMI2: Full CSDMS BMI v2.0 C surface for libswap_bmi.so.
 * All ~25 CSDMS BMI v2.0 methods implemented (no more stubs). */

/* Lifecycle */
int initialize(const char *config_file, int n);
int update(void);
int update_until(double target_time);
int finalize(void);

/* Time */
int get_current_time(double *t);
int get_start_time(double *t);
int get_end_time(double *t);
int get_time_step(double *dt);
int get_time_units(char *units_buf, int n);

/* Component info */
int get_component_name(char *name_buf, int n);
int get_input_item_count(int *count);
int get_output_item_count(int *count);
int get_input_var_names(char *buf, int n);
int get_output_var_names(char *buf, int n);

/* Variable accessors */
int get_value_double(const char *var_name, int n, double *dest);
int set_value_double(const char *var_name, int n, const double *src);

/* Variable metadata */
int get_var_type(const char *var_name, char *type_buf, int n);
int get_var_units(const char *var_name, char *units_buf, int n);
int get_var_grid(const char *var_name, int *grid_id);
int get_var_itemsize(const char *var_name, int *sz);
int get_var_nbytes(const char *var_name, int *nb);
int get_var_location(const char *var_name, char *loc_buf, int n);

/* Grid metadata — single 1D vertical soil column (grid id = 0) */
int get_grid_type(int grid_id, char *type_buf, int n);
int get_grid_rank(int grid_id, int *rank);
int get_grid_size(int grid_id, int *sz);
int get_grid_shape(int grid_id, int *shape_arr, int max_n);
int get_grid_node_count(int grid_id, int *count);
int get_grid_z(int grid_id, double *z_arr, int max_n);
