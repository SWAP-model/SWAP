/* SS-DRV Phase 1: minimal BMI C surface for libswap_bmi.so.
 * Full CSDMS BMI v2.0 compliance arrives in Phase 2. */

int initialize(const char *config_file, int n);
int update(void);
int finalize(void);
int get_current_time(double *t);
int get_time_step(double *dt);
int get_value_double(const char *var_name, int n, double *dest);
int set_value_double(const char *var_name, int n, const double *src);
