mutable struct EnergyTerms
    evalsum::Float64
    E_kin::Float64
    E_kin_core::Float64
    E_cl::Float64
    E_vcl::Float64
    E_nn::Float64
    E_en::Float64
    E_har::Float64
    E_mad::Float64
    E_vxc::Float64
    E_bxc::Float64
    E_bext::Float64
    E_xc::Float64
    E_TS::Float64
    E_tot::Float64
end

function EnergyTerms()
    evalsum = 0.0
    E_kin = 0.0
    E_kin_core = 0.0
    E_cl = 0.0
    E_vcl = 0.0
    E_nn = 0.0
    E_en = 0.0
    E_har = 0.0
    E_mad = 0.0
    E_vxc = 0.0
    E_bxc = 0.0
    E_bext = 0.0
    E_xc = 0.0
    E_TS = 0.0
    E_tot = 0.0
    return EnergyTerms(
        evalsum,
        E_kin,
        E_kin_core,
        E_cl,
        E_vcl,
        E_nn,
        E_en,
        E_har,
        E_mad,
        E_vxc,
        E_bxc,
        E_bext,
        E_xc,
        E_TS,
        E_tot
    )
end

import Base: +
function +(e1::EnergyTerms, e2::EnergyTerms)
    return EnergyTerms(
        e1.evalsum    +   e2.evalsum,
        e1.E_kin      +   e2.E_kin,
        e1.E_kin_core +   e2.E_kin_core,
        e1.E_cl       +   e2.E_cl,
        e1.E_vcl      +   e2.E_vcl,
        e1.E_nn       +   e2.E_nn,
        e1.E_en       +   e2.E_en,
        e1.E_har      +   e2.E_har,
        e1.E_mad      +   e2.E_mad,
        e1.E_vxc      +   e2.E_vxc,
        e1.E_bxc      +   e2.E_bxc,
        e1.E_bext     +   e2.E_bext,
        e1.E_xc       +   e2.E_xc,
        e1.E_TS       +   e2.E_TS,
        e1.E_tot      +   e2.E_tot
    )
end

import Base: -
function -(e1::EnergyTerms, e2::EnergyTerms)
    return EnergyTerms(
        e1.evalsum    -   e2.evalsum,
        e1.E_kin      -   e2.E_kin,
        e1.E_kin_core -   e2.E_kin_core,
        e1.E_cl       -   e2.E_cl,
        e1.E_vcl      -   e2.E_vcl,
        e1.E_nn       -   e2.E_nn,
        e1.E_en       -   e2.E_en,
        e1.E_har      -   e2.E_har,
        e1.E_mad      -   e2.E_mad,
        e1.E_vxc      -   e2.E_vxc,
        e1.E_bxc      -   e2.E_bxc,
        e1.E_bext     -   e2.E_bext,
        e1.E_xc       -   e2.E_xc,
        e1.E_TS       -   e2.E_TS,
        e1.E_tot      -   e2.E_tot
    )
end

import Base: copy
function copy(e::EnergyTerms)
    return EnergyTerms(
        e.evalsum,
        e.E_kin,
        e.E_kin_core,
        e.E_cl,
        e.E_vcl,
        e.E_nn,
        e.E_en,
        e.E_har,
        e.E_mad,
        e.E_vxc,
        e.E_bxc,
        e.E_bext,
        e.E_xc,
        e.E_TS,
        e.E_tot
    )
end

function print_info(e::EnergyTerms; prefix_str="")
    @printf("%s sum of eigenvalues      %22.12f\n", prefix_str, e.evalsum)
    @printf("%s electron kinetic        %22.12f\n", prefix_str, e.E_kin)
    @printf("%s core electron kinetic   %22.12f\n", prefix_str, e.E_kin_core)
    @printf("%s Coulomb                 %22.12f\n", prefix_str, e.E_cl)
    @printf("%s Coulomb potential       %22.12f\n", prefix_str, e.E_vcl)
    @printf("%s nuclear-nuclear         %22.12f\n", prefix_str, e.E_nn)
    @printf("%s electron-nuclear        %22.12f\n", prefix_str, e.E_en)
    @printf("%s Hartree                 %22.12f\n", prefix_str, e.E_har)
    @printf("%s Madelung                %22.12f\n", prefix_str, e.E_mad)
    @printf("%s xc potential            %22.12f\n", prefix_str, e.E_vxc)
    @printf("%s xc effective B-field    %22.12f\n", prefix_str, e.E_bxc)
    @printf("%s external B-field        %22.12f\n", prefix_str, e.E_bext)
    @printf("%s XC                      %22.12f\n", prefix_str, e.E_xc)
    @printf("%s electron entropic       %22.12f\n", prefix_str, e.E_TS)
    println("--------------------------------------------------")
    @printf("%s Total energy            %22.12f\n", prefix_str, e.E_tot)
    return
end

import Base: abs
function abs(e)
    return EnergyTerms(
        abs(e.evalsum),
        abs(e.E_kin),
        abs(e.E_kin_core),
        abs(e.E_cl),  
        abs(e.E_vcl),    
        abs(e.E_nn),  
        abs(e.E_en),   
        abs(e.E_har),   
        abs(e.E_mad),    
        abs(e.E_vxc),   
        abs(e.E_bxc),   
        abs(e.E_bext),    
        abs(e.E_xc), 
        abs(e.E_TS),   
        abs(e.E_tot)   
    )
end