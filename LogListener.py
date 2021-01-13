import logging


def log_linear_regression(soft_name, model, r_sq):
    logging.info(soft_name + f' linear regression r_2:\n{r_sq}\n')
    logging.info(soft_name + f' b0 (intercept):\n{model.intercept_}\n')
    logging.info(soft_name + f' k (slope):\n{model.coef_[0]}\n')
    logging.info(soft_name + f' factor (k^-1):\n{1.0 / model.coef_[0]}\n')
