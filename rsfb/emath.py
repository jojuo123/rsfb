"""Extended math functions for numerically stable HMM
    """

import math


class ExtMath():
    @staticmethod
    def eexp(x):
        """exteded exponential function

        Args:
            x (float): logarithm probability

        Returns:
            float: exponential of x
        """
        exp_x = 0.0
        if not math.isnan(x):
            exp_x = math.exp(x)

        return exp_x

    @staticmethod
    def eln(x):
        """extened logarithm function

        Args:
            x (float): exponential probability

        Returns:
            float: logarithm of x
        """
        ln_x = math.nan

        try:
            ln_x = math.log(x)
        except Exception:
            print("[Negative or Zero Input Error] ext_ln(x)")
            pass

        return ln_x

    @staticmethod
    def eln_sum(ln_x, ln_y):
        """extended logarithm sum function

        Args:
            x (float): exponential probability
            y (float): exponential probability

        Returns:
            float: exteded logartihm of sum of x and y
        """
        ln_sum = 0.0
        if math.isnan(ln_x) or math.isnan(ln_y):
            if math.isnan(ln_x):
                ln_sum = ln_y
            else:
                ln_sum = ln_x
        else:
            if ln_x > ln_y:
                ln_sum = ln_x + math.log(1 +
                                         math.exp(ln_y - ln_x))
            else:
                ln_sum = ln_y + math.log(1 +
                                         math.exp(ln_x - ln_y))

        return ln_sum

    @staticmethod
    def eln_product(ln_x, ln_y):
        ln_product = 0.0

        if math.isnan(ln_x) or math.isnan(ln_y):
            ln_product = math.nan
        else:
            ln_product = ln_x + ln_y

        return ln_product
