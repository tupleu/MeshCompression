#![allow(unused_imports)]
#![allow(dead_code)]
#![allow(unused_variables)]

use std::env;

mod render;
use render::{Mesh, Vertex, run};

use image::GenericImageView;

use rand;



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
	
	for i in 0..100 {
		let edge = img_mesh.get_random_edge();
		let color_diff = Mesh::color_diff(&edge);
		if color_diff == [0.0,0.0,0.0] {
			img_mesh.collapse_edge(edge);
		}
		
	}
	
	pollster::block_on(run(&img_mesh.extract_vertices(), &img_mesh.extract_indices(), wire));
}
