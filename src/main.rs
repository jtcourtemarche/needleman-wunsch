use std::fmt;
use std::cmp;

fn main() {
    let mut seq1: String = String::new();
    let mut seq2: String = String::new();

    println!("Enter two DNA sequences: ");

    let stdin = std::io::stdin();

    stdin.read_line(&mut seq1)
        .expect("Failed to read sequence #1 from buffer!");
    stdin.read_line(&mut seq2)
        .expect("Failed to read sequence #2 from buffer!");

    let mut ed = EditMatrix::new(seq1, seq2);
    println!("Before: \n{}", ed);
    let edit_distance = ed.opt_distance();
    println!("After: \n{}", ed);

    println!("Edit Distance: {}\nAlignment: \n{}", edit_distance, ed.alignment());
}

struct EditMatrix {
    first: String, 
    second: String, 

    N: usize,
    M: usize,
    matrix: Vec<Vec<i32>>,
}

impl EditMatrix {
    // construct edit distance object 
    fn new(seq1: String, seq2: String) -> Self {
        let first = seq1.trim();
        let second = seq2.trim();

        let n = first.chars().count();
        let m = second.chars().count();

        Self {
            first: String::from(first),
            second: String::from(second),
            N: n,
            M: m,
            matrix: {vec![vec![0; m + 1]; n + 1]},
        }
    }

    // returns the penalty for aligning chars a and b (returns 0 or 1)
    fn penalty(a: char, b: char) -> i32 {
        !(a == b) as i32
    } 

    // returns the min of the 3 arguments
    fn min(a: i32, b: i32, c: i32) -> i32 {
        cmp::min::<i32>(
            cmp::min::<i32>(a, b), c)
    }

    // populates the matrix based on having the two strings,
    // and returns the optimal distance
    pub fn opt_distance(&mut self) -> i32 {
        // fill out the NxM matrix per the min-of-three matrixions formula
        // bottom to top, right to left
        for i in (0..self.N).rev() {
            self.matrix[i][self.M] = ((self.N - i) * 2) as i32;
        }
        for i in (0..self.M).rev() {
            self.matrix[self.N][i] = ((self.M - i) * 2) as i32;
        }

        let mut p: i32;

        // fill the rest of the matrix using Needleman-Wunsch
        for i in (0..self.N).rev() {
            for j in (0..self.M).rev() {
                p = EditMatrix::penalty(
                    self.first.chars().nth(i).unwrap(), self.second.chars().nth(j).unwrap());

                self.matrix[i][j] = EditMatrix::min(
                    self.matrix[i+1][j+1] + p,
                    self.matrix[i+1][j] + 2,
                    self.matrix[i][j+1] + 2,
                );
            }
        }

        // initialize values on the last column & last row
        self.matrix[0][0]
    }

    // traces the matrix and returns a string that can be
    // printed to display the actual alignment (multi-line string w/ new line)
    pub fn alignment(&mut self) -> String {
        let mut result = String::new();

        let (mut i, mut j) = (0, 0);
        let mut x: Option<char> = None;
        let mut y: Option<char> = None;

        loop {
            if i >= self.N || j >= self.M {
                break;
            }

            x = Some(self.first.chars().nth(i)
                .expect("Alignment failed due to string not being able to be indexed."));
            y = Some(self.second.chars().nth(j)
                .expect("Alignment failed due to string not being able to be indexed."));

            // diagonal match 
            if self.matrix[i][j] == self.matrix[i+1][j+1] + 1 || x == y {
                result.push_str(&format!("{} {}", x.unwrap(), y.unwrap()));

                result += match x == y {
                    true => " 0",
                    false => " 1",
                };

                i += 1;
                j += 1;
            }
            // bottom match
            else if self.matrix[i][j] == self.matrix[i+1][j] + 2 {
                result.push_str(&format!("{} - 2", x.unwrap()));
                i += 1;
            }
            // right match
            else if self.matrix[i][j] == self.matrix[i][i+1] + 2 {
                result.push_str(&format!("- {} 2", y.unwrap()));
                j += 1;
            }

            result += "\n";
        }        

        // handle edge case
        if x.is_some() && y.is_some() {
            if i < self.N {
                result.push_str(&format!("{} - 2\n", x.unwrap()));
            } 
            else if j < self.M {
                result.push_str(&format!("- {} 2\n", y.unwrap()));
            }
        }

        result
    }
}

// prints out matrix
impl fmt::Display for EditMatrix {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut matrix = String::from("");
        for i in &self.matrix {
            for j in i {
                matrix = format!("{}{}\t", matrix, j);
            }
            matrix.push_str("\n");
        }
        write!(f, "{}", matrix)
    }
}

#[cfg(test)]
mod tests {
    use crate::EditMatrix;

    #[test]
    fn test_penalty() {
        assert_eq!(
            EditMatrix::penalty('a', 'a'), 0
        );
        assert_eq!(
            EditMatrix::penalty('a', 'b'), 1
        );
    }

    #[test]
    fn test_min() {
        assert_eq!(EditMatrix::min(1, 2, 3), 1);
    }

    #[test]
    fn test_opt_distance() {
        let mut test_ed = EditMatrix::new(
            String::from("AACAGTTACC"),
            String::from("TAAGGTCA")
        );
        assert_eq!(
            test_ed.opt_distance(),
            7
        )
    }

    #[test]
    fn test_alignment() {
        let mut test_ed = EditMatrix::new(
            String::from("AACAGTTACC"),
            String::from("TAAGGTCA")
        );
        test_ed.opt_distance();
        assert_eq!(
            test_ed.alignment(),
            String::from("A T 1\nA A 0\nC - 2\nA A 0\nG G 0\nT G 1\nT T 0\nA - 2\nC C 0\nC A 1\n")
        );
    }
}