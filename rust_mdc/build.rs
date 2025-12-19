use std::process::Command;
use std::env;
use std::path::Path;

fn main() {

    let target_os = std::env::var("CARGO_CFG_TARGET_OS").unwrap();
    if target_os == "windows" {
        // Assuming you have `libmdc.lib` compiled for Windows in `../bin`
        //println!("cargo:rustc-link-search=native=../x64/Release/");
        //println!("cargo:rustc-link-lib=static=mdc_lib");
        //for some reason the above two lines does not work and i need to do as follows
        println!("cargo:rustc-link-arg=../x64/Release/mdc_lib.lib");
        
        // Link against mdc and required system libraries (seems they are not needed?)
        //println!("cargo:rustc-link-lib=dylib=stdc++");
        //println!("cargo:rustc-link-lib=dylib=msvcrt");
        //println!("cargo:rustc-link-lib=dylib=ws2_32");
        //println!("cargo:rustc-link-lib=dylib=ucrt");
        //println!("cargo:rustc-link-lib=dylib=vcruntime");
        //println!("cargo:rustc-link-lib=dylib=cpp");

    }
    else {

        let repo_root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
        let lib_path = repo_root.join("bin/libmdc.a");
        
        if !lib_path.exists() {
            let _ = Command::new("make")
                .arg("-j")
                .arg("libmdc")
                .current_dir(repo_root)
                .status()
                .expect("Failed to build C++ library");
        }

        println!("cargo:rustc-link-search=native={}", repo_root.join("bin").display());
            
        // Force the linker arguments in the correct order using `cargo:rustc-link-arg`.
        // We wrap everything in --start-group ... --end-group so the linker won't discard
        // any symbols prematurely.
        println!("cargo:rustc-link-lib=static=mdc");
        //println!("cargo:rustc-link-arg=-Wl,--start-group");
        //println!("cargo:rustc-link-arg=-lmdc");
        println!("cargo:rustc-link-lib=stdc++");
        println!("cargo:rustc-link-lib=m");
        println!("cargo:rustc-link-lib=pthread");
        //println!("cargo:rustc-link-arg=-Wl,--end-group");
    }
}
