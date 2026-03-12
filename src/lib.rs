pub mod kvmer;
pub mod seeding;
pub mod types;
pub mod analyze;
pub mod calibrate;
pub mod cmdline;
pub mod inference;
pub mod utils;
pub mod sketch;
pub mod constants;
pub mod mapping;
pub mod huber;

#[cfg(target_arch = "x86_64")]
pub mod avx2_seeding;


#[cfg(test)]
mod tests {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;

    use crate::utils::*;

    #[test]
    fn test_kmer_neighbors() {
        let value = 0b11000001110111;

        let value_length = 7;
        //let key_length = 5; // arbitrary

        //let kvmer = kvmer::KVmerSet::new(key_length, value_length, false);
        //let neighbors = _get_neighbors(value, value_length, false);
        println!("Original value: {}", _kmer_to_string(value, value_length));
        _show_neighbors(value, value_length, true);
        
        assert_eq!(1., 2.);
    }
}