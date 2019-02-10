//! Common utilities shared by unit tests.
#![cfg(test)]

use std::collections::HashSet;
use std::f32;

use nalgebra::{Point3, Vector3};
use obj::raw::object::Polygon;
use obj::*;

use crate::aabb::{Bounded, AABB};
use crate::bounding_hierarchy::BoundingHierarchy;
use crate::ray::Ray;

/// A vector represented as a tuple
pub type TupleVec = (f32, f32, f32);

/// Convert a `TupleVec` to a `nalgebra` point.
pub fn tuple_to_point(tpl: &TupleVec) -> Point3<f32> {
    Point3::new(tpl.0, tpl.1, tpl.2)
}

/// Convert a `TupleVec` to a `nalgebra` vector.
pub fn tuple_to_vector(tpl: &TupleVec) -> Vector3<f32> {
    Vector3::new(tpl.0, tpl.1, tpl.2)
}

/// Define some `Bounded` structure.
pub struct UnitBox {
    pub id: i32,
    pub pos: Point3<f32>,
}

impl UnitBox {
    pub fn new(id: i32, pos: Point3<f32>) -> UnitBox {
        UnitBox { id: id, pos: pos }
    }
}

/// `UnitBox`'s `AABB`s are unit `AABB`s centered on the box's position.
impl Bounded for UnitBox {
    fn aabb(&self) -> AABB {
        let min = self.pos + Vector3::new(-0.5, -0.5, -0.5);
        let max = self.pos + Vector3::new(0.5, 0.5, 0.5);
        AABB::with_bounds(min, max)
    }
}

/// Generate 21 `UnitBox`s along the X axis centered on whole numbers (-10,9,..,10).
/// The index is set to the rounded x-coordinate of the box center.
pub fn generate_aligned_boxes() -> Vec<UnitBox> {
    // Create 21 boxes along the x-axis
    let mut shapes = Vec::new();
    for x in -10..11 {
        shapes.push(UnitBox::new(x, Point3::new(x as f32, 0.0, 0.0)));
    }
    shapes
}

/// Creates a `BoundingHierarchy` for a fixed scene structure.
pub fn build_some_bh<BH: BoundingHierarchy>() -> (Vec<UnitBox>, BH) {
    let boxes = generate_aligned_boxes();
    let bounds: Vec<AABB> = boxes.iter().map(|a| a.aabb()).collect();
    let bh = BH::build(bounds.as_slice());
    (boxes, bh)
}

/// Given a ray, a bounding hierarchy, the complete list of shapes in the scene and a list of
/// expected hits, verifies, whether the ray hits only the expected shapes.
fn traverse_and_verify<BH: BoundingHierarchy>(
    ray_origin: Point3<f32>,
    ray_direction: Vector3<f32>,
    all_shapes: &Vec<UnitBox>,
    bh: &BH,
    expected_shapes: &HashSet<i32>,
) {
    let ray = Ray::new(ray_origin, ray_direction);
    let hit_shapes = bh.traverse(&ray, all_shapes);

    assert_eq!(expected_shapes.len(), hit_shapes.len());
    for shape in hit_shapes {
        assert!(expected_shapes.contains(&shape.id));
    }
}

/// Perform some fixed intersection tests on BH structures.
pub fn traverse_some_bh<BH: BoundingHierarchy>() {
    let (all_shapes, bh) = build_some_bh::<BH>();

    {
        // Define a ray which traverses the x-axis from afar.
        let origin = Point3::new(-1000.0, 0.0, 0.0);
        let direction = Vector3::new(1.0, 0.0, 0.0);
        let mut expected_shapes = HashSet::new();

        // It should hit everything.
        for id in -10..11 {
            expected_shapes.insert(id);
        }
        traverse_and_verify(origin, direction, &all_shapes, &bh, &expected_shapes);
    }

    {
        // Define a ray which traverses the y-axis from afar.
        let origin = Point3::new(0.0, -1000.0, 0.0);
        let direction = Vector3::new(0.0, 1.0, 0.0);

        // It should hit only one box.
        let mut expected_shapes = HashSet::new();
        expected_shapes.insert(0);
        traverse_and_verify(origin, direction, &all_shapes, &bh, &expected_shapes);
    }

    {
        // Define a ray which intersects the x-axis diagonally.
        let origin = Point3::new(6.0, 0.5, 0.0);
        let direction = Vector3::new(-2.0, -1.0, 0.0);

        // It should hit exactly three boxes.
        let mut expected_shapes = HashSet::new();
        expected_shapes.insert(4);
        expected_shapes.insert(5);
        expected_shapes.insert(6);
        traverse_and_verify(origin, direction, &all_shapes, &bh, &expected_shapes);
    }
}

/// A triangle struct. Instance of a more complex `Bounded` primitive.
#[derive(Debug)]
pub struct Triangle {
    pub a: Point3<f32>,
    pub b: Point3<f32>,
    pub c: Point3<f32>,
    aabb: AABB,
    node_index: usize,
}

impl Triangle {
    pub fn new(a: Point3<f32>, b: Point3<f32>, c: Point3<f32>) -> Triangle {
        Triangle {
            a: a,
            b: b,
            c: c,
            aabb: AABB::empty().grow(&a).grow(&b).grow(&c),
            node_index: 0,
        }
    }
}

impl Bounded for Triangle {
    fn aabb(&self) -> AABB {
        self.aabb
    }
}

impl FromRawVertex for Triangle {
    fn process(
        vertices: Vec<(f32, f32, f32, f32)>,
        _: Vec<(f32, f32, f32)>,
        polygons: Vec<Polygon>,
    ) -> ObjResult<(Vec<Self>, Vec<u16>)> {
        // Convert the vertices to `Point3`s.
        let points = vertices
            .into_iter()
            .map(|v| Point3::new(v.0, v.1, v.2))
            .collect::<Vec<_>>();

        // Estimate for the number of triangles, assuming that each polygon is a triangle.
        let mut triangles = Vec::with_capacity(polygons.len());
        {
            let mut push_triangle = |indices: &Vec<usize>| {
                let mut indices_iter = indices.iter();
                let anchor = points[*indices_iter.next().unwrap()];
                let mut second = points[*indices_iter.next().unwrap()];
                for third_index in indices_iter {
                    let third = points[*third_index];
                    triangles.push(Triangle::new(anchor, second, third));
                    second = third;
                }
            };

            // Iterate over the polygons and populate the `Triangle`s vector.
            for polygon in polygons.into_iter() {
                match polygon {
                    Polygon::P(ref vec) => push_triangle(vec),
                    Polygon::PT(ref vec) | Polygon::PN(ref vec) => {
                        push_triangle(&vec.iter().map(|vertex| vertex.0).collect())
                    }
                    Polygon::PTN(ref vec) => {
                        push_triangle(&vec.iter().map(|vertex| vertex.0).collect())
                    }
                }
            }
        }
        Ok((triangles, Vec::new()))
    }
}

/// Loads the sponza model.
#[cfg(feature = "bench")]
pub fn load_sponza_scene() -> (Vec<Triangle>, AABB) {
    use std::fs::File;
    use std::io::BufReader;

    let file_input =
        BufReader::new(File::open("media/sponza.obj").expect("Failed to open .obj file."));
    let sponza_obj: Obj<Triangle> = load_obj(file_input).expect("Failed to decode .obj file data.");
    let triangles = sponza_obj.vertices;

    let mut bounds = AABB::empty();
    for triangle in &triangles {
        bounds.join_mut(&triangle.aabb());
    }

    (triangles, bounds)
}

/// Creates a `Ray` from the random `seed`. Mutates the `seed`.
/// The Ray origin will be inside the `bounds` and point to some other point inside this
/// `bounds`.
#[cfg(feature = "bench")]
pub fn create_ray(seed: &mut u64, bounds: &AABB) -> Ray {
    let origin = next_point3(seed, bounds);
    let direction = next_point3(seed, bounds).coords;
    Ray::new(origin, direction)
}

/// Benchmark the construction of a `BoundingHierarchy` with `n` triangles.
#[cfg(feature = "bench")]
fn build_n_triangles_bh<T: BoundingHierarchy>(n: usize, b: &mut ::test::Bencher) {
    let bounds = default_bounds();
    let mut triangles = create_n_cubes(n, &bounds);
    b.iter(|| {
        T::build(&mut triangles);
    });
}

/// Benchmark the construction of a `BoundingHierarchy` with 1,200 triangles.
#[cfg(feature = "bench")]
pub fn build_1200_triangles_bh<T: BoundingHierarchy>(b: &mut ::test::Bencher) {
    build_n_triangles_bh::<T>(100, b);
}

/// Benchmark the construction of a `BoundingHierarchy` with 12,000 triangles.
#[cfg(feature = "bench")]
pub fn build_12k_triangles_bh<T: BoundingHierarchy>(b: &mut ::test::Bencher) {
    build_n_triangles_bh::<T>(1_000, b);
}

/// Benchmark the construction of a `BoundingHierarchy` with 120,000 triangles.
#[cfg(feature = "bench")]
pub fn build_120k_triangles_bh<T: BoundingHierarchy>(b: &mut ::test::Bencher) {
    build_n_triangles_bh::<T>(10_000, b);
}

/// Benchmark intersecting the `triangles` list without acceleration structures.
#[cfg(feature = "bench")]
pub fn intersect_list(triangles: &[Triangle], bounds: &AABB, b: &mut ::test::Bencher) {
    let mut seed = 0;
    b.iter(|| {
        let ray = create_ray(&mut seed, &bounds);

        // Iterate over the list of triangles.
        for triangle in triangles {
            ray.intersects_triangle(&triangle.a, &triangle.b, &triangle.c);
        }
    });
}

#[cfg(feature = "bench")]
#[bench]
/// Benchmark intersecting 120,000 triangles directly.
fn bench_intersect_120k_triangles_list(b: &mut ::test::Bencher) {
    let bounds = default_bounds();
    let triangles = create_n_cubes(10_000, &bounds);
    intersect_list(&triangles, &bounds, b);
}

#[cfg(feature = "bench")]
#[bench]
/// Benchmark intersecting Sponza.
fn bench_intersect_sponza_list(b: &mut ::test::Bencher) {
    let (triangles, bounds) = load_sponza_scene();
    intersect_list(&triangles, &bounds, b);
}

/// Benchmark intersecting the `triangles` list with `AABB` checks, but without acceleration
/// structures.
#[cfg(feature = "bench")]
pub fn intersect_list_aabb(triangles: &[Triangle], bounds: &AABB, b: &mut ::test::Bencher) {
    let mut seed = 0;
    b.iter(|| {
        let ray = create_ray(&mut seed, &bounds);

        // Iterate over the list of triangles.
        for triangle in triangles {
            // First test whether the ray intersects the AABB of the triangle.
            if ray.intersects_aabb(&triangle.aabb()) {
                ray.intersects_triangle(&triangle.a, &triangle.b, &triangle.c);
            }
        }
    });
}

#[cfg(feature = "bench")]
#[bench]
/// Benchmark intersecting 120,000 triangles with preceeding `AABB` tests.
fn bench_intersect_120k_triangles_list_aabb(b: &mut ::test::Bencher) {
    let bounds = default_bounds();
    let triangles = create_n_cubes(10_000, &bounds);
    intersect_list_aabb(&triangles, &bounds, b);
}

#[cfg(feature = "bench")]
#[bench]
/// Benchmark intersecting 120,000 triangles with preceeding `AABB` tests.
fn bench_intersect_sponza_list_aabb(b: &mut ::test::Bencher) {
    let (triangles, bounds) = load_sponza_scene();
    intersect_list_aabb(&triangles, &bounds, b);
}

#[cfg(feature = "bench")]
pub fn intersect_bh<T: BoundingHierarchy>(
    bh: &T,
    triangles: &[Triangle],
    bounds: &AABB,
    b: &mut ::test::Bencher,
) {
    let mut seed = 0;
    b.iter(|| {
        let ray = create_ray(&mut seed, bounds);

        // Traverse the `BoundingHierarchy` recursively.
        let hits = bh.traverse(&ray, triangles);

        // Traverse the resulting list of positive `AABB` tests
        for triangle in &hits {
            ray.intersects_triangle(&triangle.a, &triangle.b, &triangle.c);
        }
    });
}

/// Benchmark the traversal of a `BoundingHierarchy` with `n` triangles.
#[cfg(feature = "bench")]
pub fn intersect_n_triangles<T: BoundingHierarchy>(n: usize, b: &mut ::test::Bencher) {
    let bounds = default_bounds();
    let mut triangles = create_n_cubes(n, &bounds);
    let bh = T::build(&mut triangles);
    intersect_bh(&bh, &triangles, &bounds, b)
}

/// Benchmark the traversal of a `BoundingHierarchy` with 1,200 triangles.
#[cfg(feature = "bench")]
pub fn intersect_1200_triangles_bh<T: BoundingHierarchy>(b: &mut ::test::Bencher) {
    intersect_n_triangles::<T>(100, b);
}

/// Benchmark the traversal of a `BoundingHierarchy` with 12,000 triangles.
#[cfg(feature = "bench")]
pub fn intersect_12k_triangles_bh<T: BoundingHierarchy>(b: &mut ::test::Bencher) {
    intersect_n_triangles::<T>(1_000, b);
}

/// Benchmark the traversal of a `BoundingHierarchy` with 120,000 triangles.
#[cfg(feature = "bench")]
pub fn intersect_120k_triangles_bh<T: BoundingHierarchy>(b: &mut ::test::Bencher) {
    intersect_n_triangles::<T>(10_000, b);
}
