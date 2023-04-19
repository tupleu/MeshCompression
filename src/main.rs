#![allow(unused_imports)]
#![allow(dead_code)]
#![allow(unused_variables)]

use std::env;

mod render;
use render::{Mesh, VertexPointer, run};

use image::GenericImageView;

use rand;



fn main() {
	let v = VertexPointer::new([0.0,0.0,0.0],[0.0,0.0,0.0],0);
	let v_copy = v.clone();
	v.set_color([0.0,0.0,1.0]);
	println!("{:?}",v_copy.color());
	
	
	
	
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
	let num_collapses = if args.len() > 3{
        args[3].parse::<usize>().unwrap()
    } else {
        1
    };
	let num_undos = if args.len() > 4{
        args[4].parse::<usize>().unwrap()
    } else {
        0
    };
	let img = image::open(img_path).unwrap();
	let mut img_mesh = Mesh::from_image(img);
	
	

	/*
	  
    
	
    for _ in 0..num_collapses {
        let edge = img_mesh.get_random_edge();
        let color_diff = Mesh::color_diff(&edge);
        // if color_diff == [0.0,0.0,0.0] {
            match img_mesh.collapse_edge(edge) {
                Ok(i) => continue,
                Err(e) => continue,//println!("{:?}", e),
                // Err(e) => (),
            }
        // }
    }
		*/
	
	//for _ in 0..num_undos { img_mesh.undo_nearest_edge_collapse(0.0, 0.0); }
	
    //println!("{:?}",img_mesh.tri_count());

	
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
