use core::f32;
use matrix::format::{Compressed, Conventional};
use std::{
	collections::VecDeque,
	fs::File,
	io::{BufWriter, Write},
};

// caution do not continue further, this is a mess and was written in a panic
// this includes primal/dual simplex and branch and bound with some hacky attempts
// at keeping track of variables throughout iterations

#[derive(Clone)]
pub struct VariableInfo {
	pub var_type: String,
	pub index: usize,
}

#[derive(Clone)]
struct Problem {
	var_count: i32,
	con_coef: Compressed<f32>,
	con_rhs_coef: Compressed<f32>,
	obj_coef: Compressed<f32>,
	obj_rhs: Compressed<f32>,
	section: String,
	variable_map: Vec<VariableInfo>,
}

/// Entering variable is determined by the most negative coefficient in the objective function.
pub fn get_primal_enter_var(obj_coef: &Compressed<f32>) -> i32 {
	let mut min_value = f32::INFINITY;
	let mut min_index = -1;

	for (i, &value) in obj_coef.values.iter().enumerate() {
		if value < min_value && value < 0.0 {
			min_value = value;
			min_index = i as i32;
		}
	}

	min_index
}

/// Leaving variable is determined by the minimum ratio of the right-hand side to the entering variable's coefficient.
pub fn get_primal_leaving_var(
	con_coef: &Compressed<f32>,
	rhs_coef: &Compressed<f32>,
	enter_idx: i32,
) -> i32 {
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

	min_index
}

/// Dual entering variable is determined by the lowest ratio of the objective function coefficient to the constraint coefficient.
pub fn get_dual_enter_var(obj_coef: &Compressed<f32>, leave_coef: Vec<f32>) -> i32 {
	let mut min_ratio = f32::INFINITY;
	let mut min_index = -1;

	for (i, &value) in obj_coef.values.iter().enumerate() {
		let leaving_value = leave_coef[i];
		if leaving_value >= 0.0 {
			continue;
		}

		let ratio = (value / leaving_value).abs();
		if ratio > 0.0 && ratio < min_ratio {
			min_ratio = ratio;
			min_index = i as i32;
		}
	}

	min_index
}

/// Dual leaving variable is determined by the most negative value in the right-hand side coefficients.
pub fn get_dual_leaving_var(con_rhs_coef: &Compressed<f32>) -> i32 {
	let mut min_value = f32::INFINITY;
	let mut min_index = -1;

	for (i, &value) in con_rhs_coef.values.iter().enumerate() {
		if value < 0.0 && value < min_value {
			min_value = value;
			min_index = i as i32;
		}
	}

	min_index
}

pub fn pivot_coef(
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

pub fn get_optimal_primal(
	con_coef: &mut Compressed<f32>,
	con_rhs_coef: &mut Compressed<f32>,
	obj_coef: &mut Compressed<f32>,
	rhs_coef: &mut Compressed<f32>,
) -> Result<(), ()> {
	loop {
		// maximisation checks if all coefficients in the objective function are non-negative
		if obj_coef.values.iter().all(|&x| x >= 0.0) {
			break;
		}

		let enter_idx = get_primal_enter_var(&obj_coef);
		let leaving_idx = get_primal_leaving_var(&con_coef, &con_rhs_coef, enter_idx);
		if leaving_idx == -1 {
			return Err(());
		}

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

	Ok(())
}

pub fn get_optimal_dual(
	con_coef: &mut Compressed<f32>,
	con_rhs_coef: &mut Compressed<f32>,
	obj_coef: &mut Compressed<f32>,
	rhs_coef: &mut Compressed<f32>,
	writer: &mut BufWriter<File>,
	var_map: &Vec<VariableInfo>,
) -> Result<(), ()> {
	let mut initial_pivot = false;

	loop {
		if con_rhs_coef.values.iter().all(|&v| v > -1e-6) {
			break;
		}

		let leaving_idx = get_dual_leaving_var(&con_rhs_coef);
		let filled_con_coef = Conventional::from(con_coef.clone());
		let leave_coef = (0..filled_con_coef.columns)
			.map(|j| filled_con_coef[(leaving_idx as usize, j)])
			.collect::<Vec<_>>();

		let enter_idx = get_dual_enter_var(&obj_coef, leave_coef);
		if enter_idx == -1 {
			if initial_pivot == false {
				writeln!(
					writer,
					"<> Infeasible dual problem, no entering variable found. Attempted with leaving row {}",
					leaving_idx + 1
				)
				.ok();
			}

			return Err(());
		}

		if initial_pivot == false {
			writeln!(
				writer,
				"<> Initial pivoting for dual problem with entering variable {}{} and leaving row {}",
				var_map[enter_idx as usize].var_type,
				var_map[enter_idx as usize].index + 1,
				leaving_idx + 1
			)
			.ok();
		}

		initial_pivot = true;

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

	let primal_result = get_optimal_primal(con_coef, con_rhs_coef, obj_coef, rhs_coef);
	primal_result
}

pub fn get_branch_var(
	var_count: i32,
	con_coef: &Compressed<f32>,
	con_rhs_coef: &Compressed<f32>,
) -> (i32, f32) {
	let mut min_value = f32::INFINITY;
	let mut min_index = -1;

	let basic_cols = (0..var_count as usize)
		.filter(|j| {
			let mut found_one = false;
			for i in 0..con_coef.rows {
				let coef = con_coef.get((i, *j));
				if coef == 1.0 {
					if found_one {
						return false;
					}

					found_one = true;
				} else if coef != 0.0 {
					return false;
				}
			}

			found_one
		})
		.collect::<Vec<_>>();

	for (i, &value) in con_rhs_coef.values.iter().enumerate() {
		if value.fract() == 0.0 {
			// ignore integer solutions
			continue;
		}

		for &j in &basic_cols {
			let coef = con_coef.get((i, j));
			if coef == 0.0 || coef != 1.0 {
				continue;
			}

			let frac_diff = value - value.floor();
			if frac_diff < min_value {
				min_value = frac_diff;
				min_index = j as i32;
				break;
			}
		}
	}

	(min_index, min_value)
}

pub fn branch_and_bound(
	var_count: i32,
	con_coef: &Compressed<f32>,
	con_rhs_coef: &Compressed<f32>,
	obj_coef: &Compressed<f32>,
	obj_rhs: &Compressed<f32>,
	initial_variable_map: Vec<VariableInfo>,
) -> std::io::Result<()> {
	let file = File::create("branches.txt")?;
	let mut writer = BufWriter::new(file);

	let mut queue = VecDeque::new();
	let root_problem = Problem {
		var_count,
		con_coef: con_coef.clone(),
		con_rhs_coef: con_rhs_coef.clone(),
		obj_coef: obj_coef.clone(),
		obj_rhs: obj_rhs.clone(),
		section: "0".to_string(),
		variable_map: initial_variable_map,
	};

	queue.push_back(root_problem);
	while let Some(current_problem) = queue.pop_front() {
		writeln!(
			writer,
			"=== Processing Problem {} ===",
			current_problem.section
		)?;

		let int_obj_vars = current_problem
			.con_rhs_coef
			.values
			.iter()
			.all(|&v| v.fract() == 0.0);
		if int_obj_vars {
			writeln!(
				writer,
				"Problem {}: All variables are integers, optimal solution found.",
				current_problem.section
			)?;

			writeln!(
				writer,
				"Problem {}: Objective value: {}\n",
				current_problem.section,
				current_problem.obj_rhs.get((0, 0))
			)?;
			continue;
		}

		let (branch_var_idx, branch_var_value) = get_branch_var(
			current_problem
				.variable_map
				.iter()
				.filter(|v| v.var_type == "x")
				.count() as i32,
			&current_problem.con_coef,
			&current_problem.con_rhs_coef,
		);

		if branch_var_idx == -1 {
			writeln!(
				writer,
				"Problem {}: No branching variable found, optimal solution reached.",
				current_problem.section
			)?;
			writeln!(
				writer,
				"Problem {}: Final objective value: {}\n",
				current_problem.section,
				current_problem.obj_rhs.get((0, 0))
			)?;
			continue;
		}

		writeln!(
			writer,
			"Problem {}: Branching on variable {} with value {}",
			current_problem.section,
			branch_var_idx + 1,
			branch_var_value
		)?;

		let new_var_count = current_problem.var_count + 1;

		writeln!(
			writer,
			"--- Creating branch {}.1 ({}{} ≤ {}) ---",
			current_problem.section,
			current_problem.variable_map[branch_var_idx as usize].var_type,
			current_problem.variable_map[branch_var_idx as usize].index + 1,
			branch_var_value.floor()
		)?;

		let left_problem = create_left_branch(
			&current_problem,
			branch_var_idx,
			branch_var_value,
			new_var_count,
			&mut writer,
		)?;

		if let Some(problem) = left_problem {
			queue.push_back(problem);
		}

		writeln!(
			writer,
			"--- Creating branch {}.2 ({}{} ≥ {}) ---",
			current_problem.section,
			current_problem.variable_map[branch_var_idx as usize].var_type,
			current_problem.variable_map[branch_var_idx as usize].index + 1,
			branch_var_value.ceil()
		)?;

		let right_problem = create_right_branch(
			&current_problem,
			branch_var_idx,
			branch_var_value,
			new_var_count,
			&mut writer,
		)?;

		if let Some(problem) = right_problem {
			queue.push_back(problem);
		}

		writeln!(writer)?;
	}

	writeln!(writer, "=== Branch and Bound Algorithm Completed ===")?;
	writer.flush()?;

	println!("Branch and bound completed. Results written to 'branches.txt'");
	Ok(())
}

fn create_left_branch(
	parent: &Problem,
	branch_var_idx: i32,
	branch_var_value: f32,
	new_var_count: i32,
	writer: &mut BufWriter<File>,
) -> std::io::Result<Option<Problem>> {
	let mut lt_con_coef = parent.con_coef.clone();
	lt_con_coef.resize((parent.con_coef.rows + 1, parent.con_coef.columns + 1));

	let mut lt_con_rhs_coef = parent.con_rhs_coef.clone();
	lt_con_rhs_coef.resize((parent.con_rhs_coef.rows + 1, 1));

	let mut lt_obj_coef = parent.obj_coef.clone();
	lt_obj_coef.resize((1, parent.obj_coef.columns + 1));
	let mut lt_obj_rhs = parent.obj_rhs.clone();

	let branch_var_row_idx = (0..lt_con_coef.rows)
		.find(|&i| lt_con_coef.get((i, branch_var_idx as usize)) == 1.0)
		.expect("Branching variable not found in constraints");

	// setup the new constraint
	let lt_new_con_row = (0..lt_con_coef.columns)
		.map(|j| {
			if j == branch_var_idx as usize || j == lt_con_coef.columns - 1 {
				1.0
			} else {
				0.0
			}
		})
		.collect::<Vec<_>>();

	lt_con_rhs_coef.set((lt_con_rhs_coef.rows - 1, 0), -branch_var_value);
	for (j, &value) in lt_new_con_row.iter().enumerate() {
		let branch_row_coef = lt_con_coef.get((branch_var_row_idx, j));
		lt_con_coef.set((lt_con_coef.rows - 1, j), (branch_row_coef - value) * -1.0);
	}

	let mut new_variable_map = parent.variable_map.clone();
	let x_count = new_variable_map
		.iter()
		.filter(|v| v.var_type == "x")
		.count();
	let slack_index = new_variable_map.len() - x_count;
	new_variable_map.push(VariableInfo {
		var_type: "s".to_string(),
		index: slack_index,
	});

	let result = get_optimal_dual(
		&mut lt_con_coef,
		&mut lt_con_rhs_coef,
		&mut lt_obj_coef,
		&mut lt_obj_rhs,
		writer,
		&new_variable_map,
	);

	if result.is_ok() {
		let obj_value = lt_obj_rhs.get((0, 0));
		writeln!(
			writer,
			"Problem {}.1: Found optimal solution with objective value: {:.3}",
			parent.section, obj_value
		)?;

		write!(writer, "Problem {}.1: Variable values: ", parent.section)?;
		for j in 0..lt_con_coef.columns {
			let mut one_row: Option<usize> = None;
			let mut is_basic = true;

			for i in 0..lt_con_coef.rows {
				let val = lt_con_coef.get((i, j));
				if (val - 1.0).abs() < 1e-6 {
					if one_row.is_none() {
						one_row = Some(i);
					} else {
						is_basic = false;
						break;
					}
				} else if val.abs() > 1e-6 {
					is_basic = false;
					break;
				}
			}

			let value = if is_basic {
				if let Some(row) = one_row {
					lt_con_rhs_coef.get((row, 0))
				} else {
					0.0
				}
			} else {
				0.0
			};

			let var_name = new_variable_map
				.get(j)
				.map(|v| format!("{}{}", v.var_type, v.index + 1))
				.unwrap_or_else(|| format!("var{}", j));

			write!(writer, "{} = {:.3} ", var_name, value)?;
		}

		writeln!(writer)?;

		Ok(Some(Problem {
			var_count: new_var_count,
			con_coef: lt_con_coef,
			con_rhs_coef: lt_con_rhs_coef,
			obj_coef: lt_obj_coef,
			obj_rhs: lt_obj_rhs,
			section: format!("{}.1", parent.section),
			variable_map: new_variable_map,
		}))
	} else {
		writeln!(
			writer,
			"Problem {}.1: Infeasible or unbounded",
			parent.section
		)?;
		Ok(None)
	}
}

fn create_right_branch(
	parent: &Problem,
	branch_var_idx: i32,
	branch_var_value: f32,
	new_var_count: i32,
	writer: &mut BufWriter<File>,
) -> std::io::Result<Option<Problem>> {
	let mut gt_con_coef = parent.con_coef.clone();
	gt_con_coef.resize((parent.con_coef.rows + 1, parent.con_coef.columns + 1));

	let mut gt_con_rhs_coef = parent.con_rhs_coef.clone();
	gt_con_rhs_coef.resize((parent.con_rhs_coef.rows + 1, 1));

	let mut gt_obj_coef = parent.obj_coef.clone();
	gt_obj_coef.resize((1, parent.obj_coef.columns + 1));
	let mut gt_obj_rhs = parent.obj_rhs.clone();

	let branch_var_row_idx = (0..gt_con_coef.rows)
		.find(|&i| gt_con_coef.get((i, branch_var_idx as usize)) == 1.0)
		.expect("Branching variable not found in constraints");

	// setup the new constraint
	let gt_new_con_row = (0..gt_con_coef.columns)
		.map(|j| {
			if j == branch_var_idx as usize {
				1.0
			} else if j == gt_con_coef.columns - 1 {
				-1.0
			} else {
				0.0
			}
		})
		.collect::<Vec<_>>();

	dbg!(branch_var_value);
	gt_con_rhs_coef.set((gt_con_rhs_coef.rows - 1, 0), branch_var_value - 1.0);
	for (j, &value) in gt_new_con_row.iter().enumerate() {
		let branch_row_coef = gt_con_coef.get((branch_var_row_idx, j));
		gt_con_coef.set((gt_con_coef.rows - 1, j), branch_row_coef - value);
	}

	let mut new_variable_map = parent.variable_map.clone();
	let x_count = new_variable_map
		.iter()
		.filter(|v| v.var_type == "x")
		.count();
	let excess_index = new_variable_map.len() - x_count;
	new_variable_map.push(VariableInfo {
		var_type: "e".to_string(),
		index: excess_index,
	});

	let result = get_optimal_dual(
		&mut gt_con_coef,
		&mut gt_con_rhs_coef,
		&mut gt_obj_coef,
		&mut gt_obj_rhs,
		writer,
		&new_variable_map,
	);

	if result.is_ok() {
		let obj_value = gt_obj_rhs.get((0, 0));
		writeln!(
			writer,
			"Problem {}.2: Found optimal solution with objective value: {:.3}",
			parent.section, obj_value
		)?;

		write!(writer, "Problem {}.2: Variable values: ", parent.section)?;
		for j in 0..gt_con_coef.columns {
			let mut one_row: Option<usize> = None;
			let mut is_basic = true;

			for i in 0..gt_con_coef.rows {
				let val = gt_con_coef.get((i, j));
				if (val - 1.0).abs() < 1e-6 {
					if one_row.is_none() {
						one_row = Some(i);
					} else {
						is_basic = false;
						break;
					}
				} else if val.abs() > 1e-6 {
					is_basic = false;
					break;
				}
			}

			let value = if is_basic {
				if let Some(row) = one_row {
					gt_con_rhs_coef.get((row, 0))
				} else {
					0.0
				}
			} else {
				0.0
			};

			let var_name = new_variable_map
				.get(j)
				.map(|v| format!("{}{}", v.var_type, v.index + 1))
				.unwrap_or_else(|| format!("var{}", j));

			write!(writer, "{} = {:.3} ", var_name, value)?;
		}

		writeln!(writer)?;

		Ok(Some(Problem {
			var_count: new_var_count,
			con_coef: gt_con_coef,
			con_rhs_coef: gt_con_rhs_coef,
			obj_coef: gt_obj_coef,
			obj_rhs: gt_obj_rhs,
			section: format!("{}.2", parent.section),
			variable_map: new_variable_map,
		}))
	} else {
		writeln!(
			writer,
			"Problem {}.2: Infeasible or unbounded",
			parent.section
		)?;
		Ok(None)
	}
}

pub fn print_tableau(
	con_coef: &Compressed<f32>,
	con_rhs_coef: &Compressed<f32>,
	obj_coef: &Compressed<f32>,
	rhs_coef: &Compressed<f32>,
) {
	obj_coef.values.iter().enumerate().for_each(|(i, &value)| {
		print!("|");
		print!("{:8.2} ", value);
	});

	print!("|");
	print!("{:8.2} ", rhs_coef.get((0, 0)));
	println!();

	for i in 0..con_coef.rows {
		println!();
		for j in 0..con_coef.columns {
			print!("|");
			print!("{:8.2} ", con_coef.get((i, j)));
		}

		print!("|");
		print!("{:8.2} ", con_rhs_coef.get((i, 0)));
	}

	println!();
}

pub fn print_comp_matrix(matrix: &Compressed<f32>) {
	for i in 0..matrix.rows {
		for j in 0..matrix.columns {
			print!("{:8.2} ", matrix.get((i, j)));
		}

		println!();
	}
}

pub fn print_conv_matrix(matrix: &Conventional<f32>) {
	for i in 0..matrix.rows {
		for j in 0..matrix.columns {
			print!("{:8.2} ", matrix[(i, j)]);
		}

		println!();
	}
}
