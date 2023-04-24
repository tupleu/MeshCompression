#![allow(unused_imports)]
#![allow(dead_code)]
#![allow(unused_variables)]
#![allow(unused_assignments)]
#![allow(unused_must_use)]

use std::env;

mod render;
use render::{Mesh, VertexPointer, run, RANDOM};

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
	// random seed
	let random_seed = if args.len() > 3{
        args[3].parse::<u64>().unwrap()
    } else {
        0
    };
	// collapses
	let percent_collapses = if args.len() > 4{
        args[4].parse::<f64>().unwrap()
    } else {
        0.0
    };
	// num uncollapses
	let percent_uncollapses = if args.len() > 5{
        args[5].parse::<f64>().unwrap()
    } else {
        0.0
    };
	
	let img = image::open(img_path).unwrap();
	let mut img_mesh = Mesh::from_image1(img, random_seed);

	let num_collapses = ((1.0 / 6.0) * percent_collapses * img_mesh.edge_count() as f64) as usize;
	let num_uncollapses = (percent_uncollapses * num_collapses as f64) as usize;

	img_mesh.do_n_edge_collapses(num_collapses, RANDOM);
	img_mesh.undo_n_edge_collapses_at(num_uncollapses, [0.0, 0.0, 0.0]);
	
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
