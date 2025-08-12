use knapster::{VariableInfo, branch_and_bound, get_optimal_primal, print_tableau};
use matrix::{Matrix, format::Compressed};

fn main() {
	// max z = 2x1 + 3x2 + 3x3 + 5x4 + 2x5 + 4x6
	// s.t. 11x1 + 8x2 + 6x3 + 14x4 + 10x5 + 10x6 <= 40
	// x1, x2, x3, x4, x5, x6 <= 1

	let values = vec![2.0, 3.0, 3.0, 5.0, 2.0, 4.0];
	let weights = vec![11.0, 8.0, 6.0, 14.0, 10.0, 10.0];
	let max_weight = 40.0;

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
	)
	.ok();

	println!(
		"Found initial optimal solution, with objective function value: {}",
		obj_rhs.get((0, 0))
	);
	print_tableau(&con_coef, &con_rhs_coef, &obj_coef, &obj_rhs);

	let initial_vars = vec![
		VariableInfo {
			var_type: "x".to_string(),
			index: 0,
		},
		VariableInfo {
			var_type: "x".to_string(),
			index: 1,
		},
		VariableInfo {
			var_type: "x".to_string(),
			index: 2,
		},
		VariableInfo {
			var_type: "x".to_string(),
			index: 3,
		},
		VariableInfo {
			var_type: "x".to_string(),
			index: 4,
		},
		VariableInfo {
			var_type: "x".to_string(),
			index: 5,
		},
		VariableInfo {
			var_type: "s".to_string(),
			index: 0,
		},
		VariableInfo {
			var_type: "s".to_string(),
			index: 1,
		},
		VariableInfo {
			var_type: "s".to_string(),
			index: 2,
		},
		VariableInfo {
			var_type: "s".to_string(),
			index: 3,
		},
		VariableInfo {
			var_type: "s".to_string(),
			index: 4,
		},
		VariableInfo {
			var_type: "s".to_string(),
			index: 5,
		},
		VariableInfo {
			var_type: "s".to_string(),
			index: 6,
		},
	];

	branch_and_bound(
		values.len() as i32,
		&con_coef,
		&con_rhs_coef,
		&obj_coef,
		&obj_rhs,
		initial_vars,
	)
	.ok();

	// let value_names = (0..values.len())
	// 	.map(|i| format!("x{}", i + 1))
	// 	.collect::<Vec<_>>();
}
