function run_pwscf( Ham )
    run(`rm -rfv TEMP_pwscf/\*`)
    write_pwscf( Ham, prefix_dir="TEMP_pwscf" )
    cd("./TEMP_pwscf")
    run(pipeline(`pw.x`, stdin="PWINPUT", stdout="LOG1"))
    cd("../")

    pwscf_energies = read_pwscf_etotal("TEMP_pwscf/LOG1")
    pwscf_forces = read_pwscf_forces("TEMP_pwscf/LOG1")
    return pwscf_energies, pwscf_forces
end