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
	
	let img_mesh = Mesh::from_image(img);
    println!("Triangle Count: {}", img_mesh.tri_count());
	
	// specific rules for collapsing - Brian
	// * color
	// * crease rules
	// * edge preservation
	// store collapses - Dru
	// selective refinement (paper) - Dru
	// * lookup nearest edge collapse
	// * default to center
	// screenshots/gif generation of refinement - Brian
	// 
	
	pollster::block_on(run(img_mesh, wire));
}
