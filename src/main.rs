use knapster::{get_optimal_primal, print_comp_matrix};
use matrix::{Matrix, format::Compressed};

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
		"Found initial optimal solution, with objective function value: {}",
		obj_rhs.get((0, 0))
	);
}
