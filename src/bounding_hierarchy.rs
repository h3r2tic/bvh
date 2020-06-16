//! This module defines the `BoundingHierarchy` trait.

pub use crate::aabb::Bounded as BHShape;
use crate::aabb::AABB;
use crate::ray::Ray;

/// This trait defines an acceleration structure with space partitioning.
/// This structure is used to efficiently compute ray-scene intersections.
pub trait BoundingHierarchy {
    /// Creates a new [`BoundingHierarchy`] from the `shapes` slice.
    ///
    /// # Examples
    ///
    /// ```
    /// use bvh::aabb::{AABB, Bounded};
    /// use bvh::bounding_hierarchy::BoundingHierarchy;
    /// use bvh::glam::Vec3;
    /// # use bvh::bounding_hierarchy::BHShape;
    /// # pub struct UnitBox {
    /// #     pub id: i32,
    /// #     pub pos: Vec3,
    /// #     node_index: usize,
    /// # }
    /// #
    /// # impl UnitBox {
    /// #     pub fn new(id: i32, pos: Vec3) -> UnitBox {
    /// #         UnitBox {
    /// #             id: id,
    /// #             pos: pos,
    /// #             node_index: 0,
    /// #         }
    /// #     }
    /// # }
    /// #
    /// # impl Bounded for UnitBox {
    /// #     fn aabb(&self) -> AABB {
    /// #         let min = self.pos + Vec3::new(-0.5, -0.5, -0.5);
    /// #         let max = self.pos + Vec3::new(0.5, 0.5, 0.5);
    /// #         AABB::with_bounds(min, max)
    /// #     }
    /// # }
    /// #
    /// # fn create_bhshapes() -> Vec<UnitBox> {
    /// #     let mut shapes = Vec::new();
    /// #     for i in 0..1000 {
    /// #         let position = Vec3::new(i as f32, i as f32, i as f32);
    /// #         shapes.push(UnitBox::new(i, position));
    /// #     }
    /// #     shapes
    /// # }
    ///
    /// let mut shapes = create_bhshapes();
    /// // Construct a normal `BVH`.
    /// {
    ///     use bvh::bvh::BVH;
    ///     let bvh = BVH::build(&mut shapes);
    /// }
    /// ```
    ///
    /// [`BoundingHierarchy`]: trait.BoundingHierarchy.html
    ///
    fn build(shapes: &[AABB]) -> Self;

    /// Traverses the [`BoundingHierarchy`].
    /// Returns a subset of `shapes`, in which the [`AABB`]s of the elements were hit by `ray`.
    ///
    /// # Examples
    ///
    /// ```
    /// use bvh::aabb::{AABB, Bounded};
    /// use bvh::bounding_hierarchy::BoundingHierarchy;
    /// use bvh::bvh::BVH;
    /// use bvh::glam::Vec3;
    /// use bvh::ray::Ray;
    /// # use bvh::bounding_hierarchy::BHShape;
    /// # pub struct UnitBox {
    /// #     pub id: i32,
    /// #     pub pos: Vec3,
    /// #     node_index: usize,
    /// # }
    /// #
    /// # impl UnitBox {
    /// #     pub fn new(id: i32, pos: Vec3) -> UnitBox {
    /// #         UnitBox {
    /// #             id: id,
    /// #             pos: pos,
    /// #             node_index: 0,
    /// #         }
    /// #     }
    /// # }
    /// #
    /// # impl Bounded for UnitBox {
    /// #     fn aabb(&self) -> AABB {
    /// #         let min = self.pos + Vec3::new(-0.5, -0.5, -0.5);
    /// #         let max = self.pos + Vec3::new(0.5, 0.5, 0.5);
    /// #         AABB::with_bounds(min, max)
    /// #     }
    /// # }
    /// #
    /// # fn create_bvh() -> (BVH, Vec<UnitBox>) {
    /// #     let mut shapes = Vec::new();
    /// #     for i in 0..1000 {
    /// #         let position = Vec3::new(i as f32, i as f32, i as f32);
    /// #         shapes.push(UnitBox::new(i, position));
    /// #     }
    /// #     let bvh = BVH::build(&mut shapes);
    /// #     (bvh, shapes)
    /// # }
    ///
    /// let (bvh, shapes) = create_bvh();
    ///
    /// let origin = Vec3::new(0.0, 0.0, 0.0);
    /// let direction = Vec3::new(1.0, 0.0, 0.0);
    /// let ray = Ray::new(origin, direction);
    /// let hit_shapes = bvh.traverse(&ray, &shapes);
    /// ```
    ///
    /// [`BoundingHierarchy`]: trait.BoundingHierarchy.html
    /// [`AABB`]: ../aabb/struct.AABB.html
    ///
    fn traverse<'a, Shape: BHShape>(&'a self, ray: &Ray, shapes: &'a [Shape]) -> Vec<&Shape>;

    /// Prints the [`BoundingHierarchy`] in a tree-like visualization.
    ///
    /// [`BoundingHierarchy`]: trait.BoundingHierarchy.html
    ///
    fn pretty_print(&self) {}
}
