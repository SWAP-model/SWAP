/* SS-BMI2: Python-direct C ABI for libswap_bmi.so.
 * Companion to the strict-CSDMS-BMI swap_bmi.h. */

/* Lifecycle */
int swap_set_headless(int flag);
int swap_initialize_from_toml_string(const char *buf, int n);

/* Input buffers */
int swap_attach_meteo_buffer(double *ptr, int n_days, int n_cols);

/* Array views (zero-copy) */
int swap_view_array(const char *name, double **ptr, int *n);

/* Scalar getters/setters */
int swap_get_scalar(const char *name, double *value);
int swap_set_scalar(const char *name, double value);

/* Derived summaries */
typedef struct {
    double rain;
    double evap_pot;
    double evap_act;
    double transp_pot;
    double transp_act;
    double runoff;
    double drain;
    double percolation;
    double storage_change;
    double balance_error;
} swap_water_balance_t;

int swap_get_water_balance(swap_water_balance_t *summary);
