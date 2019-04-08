#[macro_use] extern crate tera;
#[macro_use] extern crate lazy_static;

use std::collections::HashMap;

use tera::{Tera, Context, Result};
use serde_json::{Value, to_value};

use std::fs::File;
use std::io::Write;
use std::{fs, io, path::Path};
use std::io::Read;

use std::env;

lazy_static! {
    pub static ref TEMPLATES: Tera = {
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
    // read command line arguments
    let args: Vec<String> = env::args().collect();

    let template_file = &args[1];
    let json_file     = &args[2];
    let out_file      = &args[3];

    println!("template file: {}", template_file);
    println!("json file {}:"    , json_file);
    println!("out file {}:"     , out_file);

    let mut file = File::open("../../acados_template/c_templates/acados_ocp_nlp.json")?;
    let mut contents = String::new();
    file.read_to_string(&mut contents)?;
    // println!("{}", contents);
    
    // Parse the string of data into serde_json::Value.
    let v: Value = serde_json::from_str(&contents)?;


    let args: Vec<String> = env::args().collect();
    match TEMPLATES.render("acados_solver.in.c", &v) {
        Ok(s) => {
            let string_list = vec!["Foo".to_string(),"Bar".to_string()];
            let joined = string_list.join("-");
            let mut f_out = File::create("acados_solver.c").expect("Unable to create file");
            f_out.write_all(s.as_bytes())?;
        },
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
