#![allow(unused_imports)]
#![allow(dead_code)]

mod render;
use render::{Vertex, Mesh, run};

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
	let img_mesh = Mesh::from_image(10, 10);
    // println!("{:?},{:?}",img_mesh.vertices(), img_mesh.indices());
    pollster::block_on(run(img_mesh.vertices(), img_mesh.indices()));
    // let (vertices, indices) = test_vertices();
    // pollster::block_on(run(&vertices, &indices));
}
