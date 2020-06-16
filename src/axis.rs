//! Axis enum for indexing three-dimensional structures.

#![allow(unused)]
use glam::Vec3;
use std::fmt::{Display, Formatter, Result};
use std::ops::{Index, IndexMut};

struct MyType<T>(T);

/// An `Axis` in a three-dimensional coordinate system.
/// Used to access `Vec3`/`Vec3` structs via index.
///
/// # Examples
/// ```
/// use bvh::axis::Axis;
///
/// let mut position = [1.0, 0.5, 42.0];
/// position[Axis::Y] *= 4.0;
///
/// assert_eq!(position[Axis::Y], 2.0);
/// ```
///
/// `nalgebra` structures are also indexable using `Axis`.
/// For reference see [the documentation]
/// (http://nalgebra.org/doc/nalgebra/struct.Vec3.html#method.index).
///
/// ```
/// extern crate bvh;
/// extern crate glam;
///
/// use bvh::axis::Axis;
/// use glam::Vec3;
///
/// # fn main() {
/// let mut position: Vec3 = Vec3::new(1.0, 2.0, 3.0);
/// position[Axis::X] = 1000.0;
///
/// assert_eq!(position[Axis::X], 1000.0);
/// # }
/// ```
#[derive(Debug, Copy, Clone, PartialEq)]
pub enum Axis {
    /// Index of the X axis.
    X = 0,

    /// Index of the Y axis.
    Y = 1,

    /// Index of the Z axis.
    Z = 2,
}

/// Display implementation for `Axis`.
impl Display for Axis {
    fn fmt(&self, f: &mut Formatter) -> Result {
        write!(
            f,
            "{}",
            match *self {
                Axis::X => "x",
                Axis::Y => "y",
                Axis::Z => "z",
            }
        )
    }
}

/// Make slices indexable by `Axis`.
impl Index<Axis> for [f32] {
    type Output = f32;

    fn index(&self, axis: Axis) -> &f32 {
        &self[axis as usize]
    }
}

/// Make `Vec3` indexable by `Axis`.
impl Index<Axis> for Vec3 {
    type Output = f32;

    fn index(&self, axis: Axis) -> &f32 {
        let self_f32 = <Self as AsRef<[f32; 3]>>::as_ref(self);
        match axis {
            Axis::X => &self_f32[0],
            Axis::Y => &self_f32[1],
            Axis::Z => &self_f32[2],
        }
    }
}

/// Make `Vec3` indexable by `Axis`.
impl Index<Axis> for MyType<Vec3> {
    type Output = f32;

    fn index(&self, axis: Axis) -> &f32 {
        let self_f32 = <Vec3 as AsRef<[f32; 3]>>::as_ref(&self.0);
        match axis {
            Axis::X => &self_f32[0],
            Axis::Y => &self_f32[1],
            Axis::Z => &self_f32[2],
        }
    }
}

/// Make slices mutably accessible by `Axis`.
impl IndexMut<Axis> for [f32] {
    fn index_mut(&mut self, axis: Axis) -> &mut f32 {
        &mut self[axis as usize]
    }
}

/*/// Make `Vec3` mutably accessible by `Axis`.
impl IndexMut<Axis> for Vec3 {
    fn index_mut(&mut self, axis: Axis) -> &mut f32 {
        match axis {
            Axis::X => &mut self.x(),
            Axis::Y => &mut self.y(),
            Axis::Z => &mut self.z(),
        }
    }
}

/// Make `Vec3` mutably accessible by `Axis`.
impl IndexMut<Axis> for MyType<Vec3> {
    fn index_mut(&mut self, axis: Axis) -> &mut f32 {
        match axis {
            Axis::X => &mut self.0.x(),
            Axis::Y => &mut self.0.y(),
            Axis::Z => &mut self.0.z(),
        }
    }
}*/

#[cfg(test)]
mod test {
    use crate::axis::Axis;
    use quickcheck::quickcheck;

    /// Test whether accessing arrays by index is the same as accessing them by `Axis`.
    quickcheck! {
        fn test_index_by_axis(tpl: (f32, f32, f32)) -> bool {
            let a = [tpl.0, tpl.1, tpl.2];

            a[0] == a[Axis::X] && a[1] == a[Axis::Y] && a[2] == a[Axis::Z]
        }
    }

    /// Test whether arrays can be mutably set, by indexing via `Axis`.
    quickcheck! {
        fn test_set_by_axis(tpl: (f32, f32, f32)) -> bool {
            let mut a = [0.0, 0.0, 0.0];

            a[Axis::X] = tpl.0;
            a[Axis::Y] = tpl.1;
            a[Axis::Z] = tpl.2;

            a[0] == tpl.0 && a[1] == tpl.1 && a[2] == tpl.2
        }
    }
}
