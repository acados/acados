#[macro_use] extern crate tera;
#[macro_use] extern crate lazy_static;

use std::collections::HashMap;

use tera::{Tera, Context, Result};
use serde_json::{Value, to_value};

use std::fs::File;
use std::{fs, io, path::Path};
use std::io::Read;

lazy_static! {
    pub static ref TEMPLATES: Tera = {
        // let mut tera = compile_templates!("templates/*");
        let mut tera = compile_templates!("../../acados_template/c_templates_tera/*");
        tera.autoescape_on(vec![".swp"]);
        tera.register_filter("do_nothing", do_nothing_filter);
        tera
    };
}

pub fn do_nothing_filter(value: Value, _: HashMap<String, Value>) -> Result<Value> {
    let s = try_get_value!("do_nothing_filter", "value", String, value);
    Ok(to_value(&s).unwrap())
}

fn main() -> io::Result<()> {
    let mut file = File::open("../../acados_template/c_templates/acados_ocp_nlp.json")?;
    let mut contents = String::new();
    file.read_to_string(&mut contents)?;
    println!("{}", contents);
    // Parse the string of data into serde_json::Value.
    let v: Value = serde_json::from_str(&contents)?;

    // match TEMPLATES.render("hello_world.in", &v) {
    match TEMPLATES.render("acados_solver.in.c", &v) {
        Ok(s) => println!("{:?}", s),
        Err(e) => {
            println!("Error: {}", e);
            for e in e.iter().skip(1) {
                println!("Reason: {}", e);
            }
        }
    };
    println!("\n -> successfully rendered templates!");
    Ok(())
}
