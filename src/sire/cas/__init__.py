__all__ = ["Expression", "lam", "Symbol", "x", "y", "LambdaSchedule"]

from ..legacy import CAS as _CAS

from .. import use_new_api as _use_new_api

_use_new_api()

Symbol = _CAS.Symbol
Expression = _CAS.Expression

lam = Symbol("lambda")
x = Symbol("x")
y = Symbol("y")

LambdaSchedule = _CAS.LambdaSchedule


def _fix_lambdaschedule():
    if hasattr(LambdaSchedule, "__orig__get_lever_values"):
        return

    try:
        LambdaSchedule.__orig__get_lever_values = LambdaSchedule.getLeverValues
    except AttributeError:
        LambdaSchedule.__orig__get_lever_values = (
            LambdaSchedule.get_lever_values
        )

    def get_lever_values(
        obj,
        lambda_values=None,
        num_lambda=101,
        initial=1.0,
        final=2.0,
        to_pandas: bool = True,
    ):
        """
        Return the values a parameter starting at 'initial' and
        going to 'final' would take for the specified lambda values.

        Return this as a pandas DataFrame is 'to_pandas' is True

        lambda_values: list[float]
            A list of lambda values to evaluate. If this is not passed
            then the lambda values will be auto-generated

        num_lambda: int
            The number of lambda values to auto-generate if a list
            of lambda values is not passed. This defaults to 101.

        initial: float
            The initial value of the parameter (value at lambda=0).
            Defaults to 0.0

        final: float
            The final value of the parameter (value at lambda=1).
            Defaults to 1.0

        to_pandas: bool
            Whether or not to return the result as a pandas DataFrame
            (defaults to True)
        """
        if lambda_values is None:
            vals = obj.__orig__get_lever_values(
                num_lambda=num_lambda, initial=initial, final=final
            )
        else:
            vals = obj.__orig__get_lever_values(
                lambda_values=lambda_values, initial=initial, final=final
            )

        if not to_pandas:
            return vals

        import pandas as pd

        if lambda_values is None:
            vals["stage"] = obj.get_lever_stages(num_lambda=num_lambda)
        else:
            vals["stage"] = obj.get_lever_stages(lambda_values=lambda_values)

        pd_vals = {}

        pd_vals["λ"] = vals["λ"]
        pd_vals["stage"] = vals["stage"]

        for key in vals.keys():
            if key.startswith("λ") and key != "λ":
                pd_vals[key] = vals[key]

        for key in vals.keys():
            if not (key.startswith("λ") or key == "stage"):
                pd_vals[key] = vals[key]

        df = pd.DataFrame(pd_vals)
        df = df.set_index("λ")

        return df

    LambdaSchedule.get_lever_values = get_lever_values


_fix_lambdaschedule()
