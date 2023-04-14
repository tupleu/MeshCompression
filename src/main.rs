#![allow(unused_imports)]
#![allow(dead_code)]
#![allow(unused_variables)]

use std::env;

mod render;
use render::{Vertex, Mesh, TriangleIndex, run};

use image::GenericImageView;

use rand;


fn test_vertices() -> (Vec<Vertex>, Vec<u16>) {
    let vertices: &[Vertex] = &[
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

    let indices: &[u16] = &[
        0, 1, 5,
        1, 2, 3,
        3, 4, 5,
        1, 3, 6,
        3, 5, 6,
        5, 1, 6,
    ];
    
    (vertices.to_vec(), indices.to_vec())
}
fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        panic!("no image path supplied")
    }
    let img_path = &args[1];
    let wire = if args.len() > 2{
        "wire".eq(&args[2])
    } else {
        false
    };
	let img = image::open(img_path).unwrap();
	// let _img1 = image::open("./tests/mario.webp").unwrap();
	// let _img1 = image::open("./tests/mariobros.png").unwrap();

	//let _img2 = image::open("./../../tests/test2-small.jpeg").unwrap();
	//let _img3 = image::open("./../../tests/test3.jpg").unwrap();
	//let _img4 = image::open("./../../tests/test4.jpg").unwrap();
	
	
	
	
	let mut img_mesh = Mesh::from_image(img);
    
	
	let mut remove_idx = vec![];
	let mut subdivide_idx = vec![];
	for tri in img_mesh.triangles() {
		let r: f32 = rand::random();
		if r > 0.9 {
			subdivide_idx.push(tri.index());
		}
	}
	for idx in subdivide_idx {
		img_mesh.subdivide_triangle(idx);
	}
	for tri in img_mesh.triangles() {
		let r: f32 = rand::random();
		if r > 0.9 {
			remove_idx.push(tri.index());
		}
	}
	for idx in remove_idx {
		img_mesh.remove_triangle(idx);
	}
	
	
	pollster::block_on(run(img_mesh.vertices(), &img_mesh.indices(), wire));
}
