//! Matrix library.

#![no_std]

extern crate alloc;

use alloc::boxed::Box;

use crate::dense::RowMajorMatrix;

pub mod dense;
pub mod mul;
pub mod sparse;
pub mod stack;

pub trait Matrix<T> {
    fn width(&self) -> usize;

    fn height(&self) -> usize;
}

/// A `Matrix` that supports randomly accessing particular coefficients.
pub trait MatrixGet<T> {
    fn get(&self, r: usize, c: usize) -> T;
}

/// A `Matrix` that supports randomly accessing particular rows.
pub trait MatrixRows<'a, T: 'a>: Matrix<T> {
    type Row: IntoIterator<Item = &'a T>;

    fn row(&'a self, r: usize) -> Self::Row;

    fn first_row(&'a self) -> Self::Row {
        self.row(0)
    }

    fn last_row(&'a self) -> Self::Row {
        self.row(self.height() - 1)
    }

    fn to_row_major_matrix(self) -> RowMajorMatrix<T>
    where
        Self: Sized,
        T: Clone,
    {
        todo!()
    }
}

impl<T> Matrix<T> for Box<dyn Matrix<T>> {
    fn width(&self) -> usize {
        self.as_ref().width()
    }

    fn height(&self) -> usize {
        self.as_ref().height()
    }
}
