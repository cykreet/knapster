use core::f32;

use matrix::{
	Matrix,
	format::{Compressed, Conventional},
};

fn get_enter_var(obj_coef: &Compressed<f32>) -> (i32, f32) {
	let mut min_value = f32::INFINITY;
	let mut min_index = -1;

	for (i, &value) in obj_coef.values.iter().enumerate() {
		if value < min_value && value < 0.0 {
			min_value = value;
			min_index = i as i32;
		}
	}

	(min_index, min_value)
}

fn get_leaving_var(
	con_coef: &Compressed<f32>,
	rhs_coef: &Compressed<f32>,
	enter_idx: i32,
) -> (i32, f32) {
	let mut min_ratio = f32::INFINITY;
	let mut min_index = -1;

	for (i, &rhs) in rhs_coef.values.iter().enumerate() {
		if rhs > 0.0 {
			let row_coef = con_coef.get((i, enter_idx as usize));
			let ratio = rhs / row_coef;
			if ratio < min_ratio {
				if ratio == 0.0 && row_coef != 1.0 {
					continue;
				}

				min_ratio = ratio;
				min_index = i as i32;
			}
		}
	}

	(min_index, min_ratio)
}

fn pivot_coef(
	con_coef: &mut Compressed<f32>,
	con_rhs_coef: &mut Compressed<f32>,
	obj_coef: &mut Compressed<f32>,
	rhs_coef: &mut Compressed<f32>,
	enter_idx: i32,
	leaving_idx: i32,
) {
	let rhs_value = con_rhs_coef.get((leaving_idx as usize, 0));
	let filled_con_coef = Conventional::from(con_coef.clone());

	// get row values for the leaving variable
	let pivot_row = (0..filled_con_coef.columns)
		.map(|j| filled_con_coef[(leaving_idx as usize, j)])
		.collect::<Vec<_>>();
	let pivot_value = pivot_row[enter_idx as usize];

	for j in 0..con_coef.columns {
		if j != enter_idx as usize {
			let new_value = pivot_row[j] / pivot_value;
			con_coef.set((leaving_idx as usize, j), new_value);
		}
	}

	con_coef.set((leaving_idx as usize, enter_idx as usize), 1.0);
	// set corresponding rhs value in pivot row
	con_rhs_coef.set((leaving_idx as usize, 0), rhs_value / pivot_value);

	for i in 0..con_coef.rows {
		if i == leaving_idx as usize {
			continue;
		}

		let factor = con_coef.get((i, enter_idx as usize));
		for j in 0..con_coef.columns {
			let pivot_coef_val = con_coef.get((leaving_idx as usize, j));
			con_coef.set((i, j), con_coef.get((i, j)) - factor * pivot_coef_val);
		}

		let pivot_rhs_val = con_rhs_coef.get((leaving_idx as usize, 0));
		con_rhs_coef.set((i, 0), con_rhs_coef.get((i, 0)) - factor * pivot_rhs_val);
	}

	let obj_factor = obj_coef.get((0, enter_idx as usize));
	for j in 0..con_coef.columns {
		let pivot_coef_val = con_coef.get((leaving_idx as usize, j));
		obj_coef.set((0, j), obj_coef.get((0, j)) - obj_factor * pivot_coef_val);
	}

	let pivot_rhs_val = con_rhs_coef.get((leaving_idx as usize, 0));
	rhs_coef.set((0, 0), rhs_coef.get((0, 0)) - obj_factor * pivot_rhs_val)
}

fn get_optimal_primal(
	con_coef: &mut Compressed<f32>,
	con_rhs_coef: &mut Compressed<f32>,
	obj_coef: &mut Compressed<f32>,
	rhs_coef: &mut Compressed<f32>,
) {
	loop {
		// maximisation checks if all coefficients in the objective function are non-negative
		if obj_coef.values.iter().all(|&x| x >= 0.0) {
			println!("Optimal solution found, breaking.");
			break;
		}

		let (enter_idx, enter_var) = get_enter_var(&obj_coef);
		let (leaving_idx, leaving_var) = get_leaving_var(&con_coef, &con_rhs_coef, enter_idx);

		// pivot the tableau
		pivot_coef(
			con_coef,
			con_rhs_coef,
			obj_coef,
			rhs_coef,
			enter_idx,
			leaving_idx,
		);
	}
}

fn print_comp_matrix(matrix: &Compressed<f32>) {
	for i in 0..matrix.rows {
		for j in 0..matrix.columns {
			print!("{:8.2} ", matrix.get((i, j)));
		}
		println!();
	}
}

fn print_conv_matrix(matrix: &Conventional<f32>) {
	for i in 0..matrix.rows {
		for j in 0..matrix.columns {
			print!("{:8.2} ", matrix[(i, j)]);
		}
		println!();
	}
}

fn main() {
	let values = vec![2.0, 3.0, 3.0, 5.0, 2.0, 4.0];
	let weights = vec![11.0, 8.0, 6.0, 14.0, 10.0, 10.0];
	let max_weight = 40.0;

	// max z = 2x1 + 2x2 + 3x3 + 4x4 + 5x5
	// s.t. 11x1 + 12x2 + 13x3 + 14x4 + 15x5 <= 40
	// x1, x2, x3, x4, x5 <= 1

	let mut obj_coef = Compressed::<f32>::zero((1, values.len() * 2 + 1));
	let mut obj_rhs = Compressed::<f32>::zero((1, 1));
	// row for each constraint, in this case each the total weight constraint, and each variable has
	// an associated <= 1 constraint
	let mut con_coef = Compressed::<f32>::zero((values.len() + 1, values.len() * 2 + 1));
	let mut con_rhs_coef = Compressed::<f32>::zero((values.len() + 1, 1));

	// coefficients for the total weight constraint
	con_rhs_coef.set((0, 0), max_weight);

	(0..values.len()).for_each(|i| {
		obj_coef.set((0, i), -values[i]);
		// rhs coefficients for the <= 1 constraints
		con_rhs_coef.set((i + 1, 0), 1.0);

		con_coef.set((0, i), weights[i]);
		// constraint variable coefficients for the <= 1 constraints
		con_coef.set((i + 1, i), 1.0);
	});

	(0..values.len() + 1).for_each(|i| {
		// slack variable coefficients
		con_coef.set((i, values.len() + i), 1.0);
	});

	get_optimal_primal(
		&mut con_coef,
		&mut con_rhs_coef,
		&mut obj_coef,
		&mut obj_rhs,
	);

	// let value_names = (0..values.len())
	// 	.map(|i| format!("x{}", i + 1))
	// 	.collect::<Vec<_>>();

	print_comp_matrix(&con_coef);
	println!(
		"Found optimal solution, with objective function value: {}",
		obj_rhs.get((0, 0))
	);
}
