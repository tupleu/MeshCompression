#![allow(unused_imports)]
#![allow(dead_code)]

mod render;
use render::{Vertex, Mesh, run};

use image::GenericImageView;

const VERTICES: &[Vertex] = &[
    Vertex {
        position: [-0.5, -0.5, 0.0],
        // color: [0.0, 0.0, 0.0],
        color: [0.0, 1.0, 0.0],
    },
    Vertex {
        position: [0.0, -0.5, 0.0],
        color: [1.0, 0.0, 0.0],
    },
    Vertex {
        position: [0.5, -0.5, 0.0],
        // color: [0.0, 0.0, 0.0],
        color: [0.0, 0.0, 1.0],
    },
    Vertex {
        position: [0.25, 0.0, 0.0],
        color: [0.0, 1.0, 0.0],
    },
    Vertex {
        position: [0.0, 0.5, 0.0],
        // color: [0.0, 0.0, 0.0],
        color: [1.0, 0.0, 0.0],
    },
    Vertex {
        position: [-0.25, 0.0, 0.0],
        color: [0.0, 0.0, 1.0],
    },
    Vertex {
        position: [0.0, -0.225, 0.0],
        // color: [1.0, 1.0, 1.0],
        color: [0.0, 0.0, 0.0],
    },
];

const INDICES: &[u16] = &[
    0, 1, 5,
    1, 2, 3,
    3, 4, 5,
    1, 3, 6,
    3, 5, 6,
    5, 1, 6,
];

fn main() {
	let _img1 = image::open("./../../tests/test1.jpg").unwrap();
	//let _img2 = image::open("./../../tests/test2-small.jpeg").unwrap();
	//let _img3 = image::open("./../../tests/test3.jpg").unwrap();
	//let _img4 = image::open("./../../tests/test4.jpg").unwrap();
	
	let img_mesh = Mesh::from_image(_img1);
    //pollster::block_on(run(img_mesh.verticies(), img_mesh.indicies()));
    pollster::block_on(run(VERTICES, INDICES));
}
