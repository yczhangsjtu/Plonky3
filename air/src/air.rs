use core::ops::{Add, Mul, Sub};

use p3_field::{AbstractExtensionField, AbstractField, AbstractionOf, ExtensionField, Field};
use p3_matrix::dense::RowMajorMatrix;
use p3_matrix::MatrixRows;

pub trait Air<AB: AirBuilder>: Sync {
    fn eval(&self, builder: &mut AB);

    fn preprocessed_trace(&self) -> Option<RowMajorMatrix<AB::F>> {
        None
    }
}

pub trait AirBuilder: Sized {
    type F: Field;

    type Expr: AbstractionOf<Self::F>
        + Add<Self::Var, Output = Self::Expr>
        + Sub<Self::Var, Output = Self::Expr>
        + Mul<Self::Var, Output = Self::Expr>;

    type Var: Into<Self::Expr>
        + Copy
        + Add<Self::F, Output = Self::Expr>
        + Add<Self::Var, Output = Self::Expr>
        + Add<Self::Expr, Output = Self::Expr>
        + Sub<Self::F, Output = Self::Expr>
        + Sub<Self::Var, Output = Self::Expr>
        + Sub<Self::Expr, Output = Self::Expr>
        + Mul<Self::F, Output = Self::Expr>
        + Mul<Self::Var, Output = Self::Expr>
        + Mul<Self::Expr, Output = Self::Expr>;

    type M: for<'a> MatrixRows<'a, Self::Var, Row = &'a [Self::Var]>;

    fn main(&self) -> Self::M;

    fn is_first_row(&self) -> Self::Expr;
    fn is_last_row(&self) -> Self::Expr;
    fn is_transition(&self) -> Self::Expr {
        self.is_transition_window(2)
    }
    fn is_transition_window(&self, size: usize) -> Self::Expr;

    /// Returns a sub-builder whose constraints are enforced only when `condition` is nonzero.
    fn when<I: Into<Self::Expr>>(&mut self, condition: I) -> FilteredAirBuilder<Self> {
        FilteredAirBuilder {
            inner: self,
            condition: condition.into(),
        }
    }

    /// Returns a sub-builder whose constraints are enforced only when `x != y`.
    fn when_ne<I1: Into<Self::Expr>, I2: Into<Self::Expr>>(
        &mut self,
        x: I1,
        y: I2,
    ) -> FilteredAirBuilder<Self> {
        self.when(x.into() - y.into())
    }

    /// Returns a sub-builder whose constraints are enforced only on the first row.
    fn when_first_row(&mut self) -> FilteredAirBuilder<Self> {
        self.when(self.is_first_row())
    }

    /// Returns a sub-builder whose constraints are enforced only on the last row.
    fn when_last_row(&mut self) -> FilteredAirBuilder<Self> {
        self.when(self.is_last_row())
    }

    /// Returns a sub-builder whose constraints are enforced on all rows except the last.
    fn when_transition(&mut self) -> FilteredAirBuilder<Self> {
        self.when(self.is_transition())
    }

    /// Returns a sub-builder whose constraints are enforced on all rows except the last `size - 1`.
    fn when_transition_window(&mut self, size: usize) -> FilteredAirBuilder<Self> {
        self.when(self.is_transition_window(size))
    }

    fn assert_zero<I: Into<Self::Expr>>(&mut self, x: I);

    fn assert_one<I: Into<Self::Expr>>(&mut self, x: I) {
        self.assert_zero(x.into() - Self::Expr::ONE);
    }

    fn assert_eq<I1: Into<Self::Expr>, I2: Into<Self::Expr>>(&mut self, x: I1, y: I2) {
        self.assert_zero(x.into() - y.into());
    }

    /// Assert that `x` is a boolean, i.e. either 0 or 1.
    fn assert_bool<I: Into<Self::Expr>>(&mut self, x: I) {
        let x = x.into();
        self.assert_zero(x.clone() * (x - Self::Expr::ONE));
    }

    fn assert_zero_ext<ExprExt, I>(&mut self, x: I)
    where
        ExprExt: AbstractExtensionField<Self::Expr>,
        I: Into<ExprExt>,
    {
        for xb in x.into().as_base_slice().iter().cloned() {
            self.assert_zero(xb);
        }
    }

    fn assert_eq_ext<ExprExt, I1, I2>(&mut self, x: I1, y: I2)
    where
        ExprExt: AbstractExtensionField<Self::Expr>,
        I1: Into<ExprExt>,
        I2: Into<ExprExt>,
    {
        self.assert_zero_ext::<ExprExt, ExprExt>(x.into() - y.into());
    }

    fn assert_one_ext<ExprExt, I>(&mut self, x: I)
    where
        ExprExt: AbstractExtensionField<Self::Expr>,
        I: Into<ExprExt>,
    {
        let xe: ExprExt = x.into();
        let parts = xe.as_base_slice();
        self.assert_one(parts[0].clone());
        for part in &parts[1..] {
            self.assert_zero(part.clone());
        }
    }
}

pub trait PairBuilder: AirBuilder {
    fn preprocessed(&self) -> Self::M;
}

pub trait PermutationAirBuilder: AirBuilder {
    type EF: ExtensionField<Self::F>;

    type ExprEF: AbstractionOf<Self::EF>
        + AbstractExtensionField<Self::Expr>
        + Add<Self::VarEF, Output = Self::ExprEF>
        + Sub<Self::VarEF, Output = Self::ExprEF>
        + Mul<Self::VarEF, Output = Self::ExprEF>;

    type VarEF: Into<Self::ExprEF>
        + Copy
        + Add<Self::EF, Output = Self::ExprEF>
        + Add<Self::VarEF, Output = Self::ExprEF>
        + Add<Self::ExprEF, Output = Self::ExprEF>
        + Sub<Self::EF, Output = Self::ExprEF>
        + Sub<Self::VarEF, Output = Self::ExprEF>
        + Sub<Self::ExprEF, Output = Self::ExprEF>
        + Mul<Self::EF, Output = Self::ExprEF>
        + Mul<Self::VarEF, Output = Self::ExprEF>
        + Mul<Self::ExprEF, Output = Self::ExprEF>;

    type MP: for<'a> MatrixRows<'a, Self::VarEF, Row = &'a [Self::VarEF]>;

    fn permutation(&self) -> Self::MP;

    fn permutation_randomness(&self) -> &[Self::EF];
}

pub struct FilteredAirBuilder<'a, AB: AirBuilder> {
    inner: &'a mut AB,
    condition: AB::Expr,
}

impl<'a, AB: PermutationAirBuilder> PermutationAirBuilder for FilteredAirBuilder<'a, AB> {
    type EF = AB::EF;
    type VarEF = AB::VarEF;
    type ExprEF = AB::ExprEF;
    type MP = AB::MP;

    fn permutation(&self) -> Self::MP {
        self.inner.permutation()
    }

    fn permutation_randomness(&self) -> &[Self::EF] {
        self.inner.permutation_randomness()
    }
}

impl<'a, AB: AirBuilder> AirBuilder for FilteredAirBuilder<'a, AB> {
    type F = AB::F;
    type Expr = AB::Expr;
    type Var = AB::Var;
    type M = AB::M;

    fn main(&self) -> Self::M {
        self.inner.main()
    }

    fn is_first_row(&self) -> Self::Expr {
        self.inner.is_first_row()
    }

    fn is_last_row(&self) -> Self::Expr {
        self.inner.is_last_row()
    }

    fn is_transition_window(&self, size: usize) -> Self::Expr {
        self.inner.is_transition_window(size)
    }

    fn assert_zero<I: Into<Self::Expr>>(&mut self, x: I) {
        self.inner.assert_zero(self.condition.clone() * x.into());
    }
}

#[cfg(test)]
mod tests {
    use p3_matrix::MatrixRows;

    use crate::{Air, AirBuilder};

    struct FibonacciAir;

    impl<AB: AirBuilder> Air<AB> for FibonacciAir {
        fn eval(&self, builder: &mut AB) {
            let main = builder.main();

            let x_0 = main.row(0)[0];
            let x_1 = main.row(1)[0];
            let x_2 = main.row(2)[0];

            builder.when_first_row().assert_zero(x_0);
            builder.when_first_row().assert_one(x_1);
            builder.when_transition().assert_eq(x_0 + x_1, x_2);
        }
    }
}
