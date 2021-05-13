/*
    Edit distance 
    by Jacob Courtemarche
*/

use std::fmt;
use std::cmp;

fn main() {
    let mut seq1: String = String::new();
    let mut seq2: String = String::new();

    let stdin = std::io::stdin();

    stdin.read_line(&mut seq1)
        .expect("Failed to read sequence #1 from buffer!");
    stdin.read_line(&mut seq2)
        .expect("Failed to read sequence #2 from buffer!");

    let mut ed = EditDistance::new(seq1, seq2);
    println!("Before: \n{}", ed);
    let edit_distance = ed.opt_distance();
    println!("After: \n{}", ed);

    println!("Edit Distance: {}\nAlignment: \n{}", edit_distance, ed.alignment());
}

struct EditDistance {
    first: String,
    second: String, 

    N: usize,
    M: usize,
    opt: Vec<Vec<i32>>,
}

impl EditDistance {
    // construct edit distance object 
    fn new(seq1: String, seq2: String) -> EditDistance {
        let first = seq1.trim();
        let second = seq2.trim();

        let n = first.chars().count();
        let m = second.chars().count();
        /*
        println!("{}{}", seq1, seq2);
        println!("{}, {}", n, m);
        */

        EditDistance{
            first: String::from(first),
            second: String::from(second),
            N: n,
            M: m,
            opt: {vec![vec![0; m + 1]; n + 1]},
        }
    }

    // returns the penalty for aligning chars a and b (returns 0 or 1)
    fn penalty(a: char, b: char) -> i32 {
        match a == b {
            true => 0,
            false => 1,
        }
    } 

    // returns the min of the 3 arguments
    fn min(a: i32, b: i32, c: i32) -> i32 {
        cmp::min::<i32>(
            cmp::min::<i32>(a, b), c)
    }

    // populates the matrix based on having the two strings,
    // and returns the optimal distance
    pub fn opt_distance(&mut self) -> i32 {
        // fill out the NxM matrix per the min-of-three options formula
        // bottom to top, right to left
        for i in (0..self.N).rev() {
            self.opt[i][self.M] = ((self.N - i) * 2) as i32;
        }
        for i in (0..self.M).rev() {
            self.opt[self.N][i] = ((self.M - i) * 2) as i32;
        }

        let mut p: i32;

        // fill the rest of the matrix using Needleman-Wunsch
        for i in (0..self.N).rev() {
            for j in (0..self.M).rev() {
                p = EditDistance::penalty(
                    self.first.chars().nth(i).unwrap(), self.second.chars().nth(j).unwrap());

                self.opt[i][j] = EditDistance::min(
                    self.opt[i+1][j+1] + p,
                    self.opt[i+1][j] + 2,
                    self.opt[i][j+1] + 2,
                );
            }
        }

        // initialize values on the last column & last row
        self.opt[0][0]
    }

    // traces the matrix and returns a string that can be
    // printed to display the actual alignment (multi-line string w/ new line)
    pub fn alignment(&mut self) -> String {
        let mut result = String::new();

        let (mut i, mut j) = (0, 0);
        let (mut x, mut y);

        loop {
            if i >= self.N || j >= self.M {
                break;
            }

            x = self.first.chars().nth(i)
                .expect("Alignment failed due to string not being able to be indexed.");
            y = self.second.chars().nth(j)
                .expect("Alignment failed due to string not being able to be indexed.");

            // diagonal match 
            if self.opt[i][j] == self.opt[i+1][j+1] + 1 || x == y {
                result.push_str(&format!("{} {}", x, y));

                result += match x == y {
                    true => " 0",
                    false => " 1",
                };

                i += 1;
                j += 1;
            }
            // bottom match
            else if self.opt[i][j] == self.opt[i+1][j] + 2 {
                result.push_str(&format!("{} - 2", x));
                i += 1;
            }
            // right match
            else if self.opt[i][j] == self.opt[i][i+1] + 2 {
                result.push_str(&format!("- {} 2", y));
                j += 1;
            }

            // handle edge case


            result += "\n";
        }        

        result
    }
}

// prints out matrix
impl fmt::Display for EditDistance {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut matrix = String::from("");
        for i in &self.opt {
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
    use crate::EditDistance;

    #[test]
    fn test_penalty() {
        assert_eq!(
            EditDistance::penalty('a', 'a'), 0
        );
        assert_eq!(
            EditDistance::penalty('a', 'b'), 1
        );
    }
}