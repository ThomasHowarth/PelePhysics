"""Production functions."""

from math import isclose, log, pi, sqrt

import cantera as ct
import symengine as sme

import ceptr.constants as cc
import ceptr.utilities as cu
import ceptr.writer as cw


def production_rate(
    fstream,
    mechanism,
    species_info,
    reaction_info,
    syms=None,
):
    """Write production rate."""
    n_species = species_info.n_species
    n_qss_species = species_info.n_qssa_species
    n_reactions = mechanism.n_reactions

    assert len(reaction_info.index) == 7

    itroe = reaction_info.index[0:2]
    isri = reaction_info.index[1:3]
    ilindemann = reaction_info.index[2:4]
    # i3body = reaction_info.index[3:5]
    # isimple = reaction_info.index[4:6]
    # ispecial = reaction_info.index[5:7]

    ntroe = itroe[1] - itroe[0]
    nsri = isri[1] - isri[0]
    nlindemann = ilindemann[1] - ilindemann[0]
    # n3body = i3body[1] - i3body[0]
    # nsimple = isimple[1] - isimple[0]
    # nspecial = ispecial[1] - ispecial[0]

    # qdot
    cw.writer(fstream)
    cw.writer(
        fstream,
        cw.comment(
            "GPU version of productionRate: no more use of thermo namespace vectors "
        ),
    )
    if n_reactions > 0:
        if n_qss_species > 0:
            cw.writer(
                fstream,
                "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void"
                " comp_qfqr(amrex::Real *  qf, amrex::Real * qr, const"
                " amrex::Real * sc, const amrex::Real * sc_qss, const"
                " amrex::Real T, const amrex::Real invT, const amrex::Real logT)",
            )
        else:
            cw.writer(
                fstream,
                "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void"
                " comp_qfqr(amrex::Real *  qf, amrex::Real * qr, const"
                " amrex::Real * sc, const amrex::Real * /*sc_qss*/,const"
                " amrex::Real T, const amrex::Real invT, const amrex::Real logT)",
            )
    else:
        cw.writer(
            fstream,
            "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void comp_qfqr(amrex::Real *"
            " /*qf*/, amrex::Real * /*qr*/, const amrex::Real * /*sc*/, const"
            " amrex::Real * /*sc_qss*/, const amrex::Real /*T*/, const amrex::Real"
            " /*invT*/, const amrex::Real /*logT*/)",
        )
    cw.writer(fstream, "{")

    if n_reactions > 0:
        # nclassd = n_reactions - nspecial
        # nCorr   = n3body + ntroe + nsri + nlindemann

        # reacs are sorted here
        for orig_idx, idx in reaction_info.idxmap.items():
            cw.writer(fstream)
            reaction = mechanism.reaction(orig_idx)
            cw.writer(
                fstream,
                cw.comment(f"reaction {orig_idx}: {reaction.equation}"),
            )
            if bool(reaction.orders):
                string = cu.qss_sorted_phase_space(
                    mechanism, species_info, reaction, reaction.orders
                )
                cw.writer(
                    fstream,
                    f"qf[{idx}] = {string};",
                )
            else:
                string = cu.qss_sorted_phase_space(
                    mechanism,
                    species_info,
                    reaction,
                    reaction.reactants,
                )
                cw.writer(
                    fstream,
                    f"qf[{idx}] = {string};",
                )

            if reaction.reversible:
                string = cu.qss_sorted_phase_space(
                    mechanism,
                    species_info,
                    reaction,
                    reaction.products,
                )
                cw.writer(
                    fstream,
                    f"qr[{idx}] = {string};",
                )
            else:
                cw.writer(fstream, f"qr[{idx}] = 0.0;")

        cw.writer(fstream)

        # Mixt concentration for PD & TB
        cw.writer(fstream, cw.comment("compute the mixture concentration"))
        cw.writer(fstream, "amrex::Real mixture = 0.0;")
        cw.writer(fstream, f"for (int i = 0; i < {n_species}; ++i) {{")
        cw.writer(fstream, "mixture += sc[i];")
        cw.writer(fstream, "}")
        cw.writer(fstream)
        if species_info.n_qssa_species > 0:
            cw.writer(
                fstream,
                f"for (int i = 0; i < {species_info.n_qssa_species}; ++i)",
            )
            cw.writer(fstream, "{")
            cw.writer(fstream, "mixture += sc_qss[i];")
            cw.writer(fstream, "}")
            cw.writer(fstream)

        # Kc stuff
        cw.writer(fstream, cw.comment("compute the Gibbs free energy"))
        cw.writer(fstream, f"amrex::Real g_RT[{species_info.n_species}];")
        cw.writer(fstream, "gibbs(g_RT, T);")
        if species_info.n_qssa_species > 0:
            cw.writer(
                fstream,
                f"amrex::Real g_RT_qss[{species_info.n_qssa_species}];",
            )
            cw.writer(fstream, "gibbs_qss(g_RT_qss, T);")
        cu.write_gibbs_exp(fstream, mechanism, species_info)

        cw.writer(fstream)

        cw.writer(
            fstream,
            cw.comment("reference concentration: P_atm / (RT) in inverse mol/m^3"),
        )
        cw.writer(
            fstream,
            f"amrex::Real refC = {cc.Patm_pa:g} /"
            f" {cc.R.to(cc.ureg.joule / (cc.ureg.mole / cc.ureg.kelvin)).m:g} *"
            " invT;",
        )
        cw.writer(fstream, "amrex::Real refCinv = 1 / refC;")

        cw.writer(fstream)

        # kfs
        cw.writer(fstream, cw.comment("Evaluate the kfs"))
        cw.writer(fstream, "amrex::Real k_f, Corr;")
        if ntroe > 0:
            cw.writer(
                fstream,
                "amrex::Real redP, F, logPred, logFcent, troe_c, troe_n, troe, F_troe;",
            )
        if nsri > 0:
            cw.writer(fstream, "amrex::Real redP, F, X, F_sri;")
        cw.writer(fstream)

        # Loop like you're going through them in the mech.Linp order
        for orig_idx, idx in sorted(reaction_info.idxmap.items()):
            reaction = mechanism.reaction(orig_idx)

            kc_e = cu.kc_exp(mechanism, species_info, reaction, -1)
            kc_conv_inv = cu.fkc_conv_inv(mechanism, species_info, reaction)

            if bool(reaction.orders):
                dim = cu.phase_space_units(reaction.orders)
            else:
                dim = cu.phase_space_units(reaction.reactants)
            third_body = reaction.third_body is not None
            if reaction.reaction_type == "pressure-dependent-Arrhenius":
                plog = True
            elif reaction.reaction_type == "three-body-pressure-dependent-Arrhenius":
                if not third_body:
                    raise ValueError(
                        f"Inconsistent reaction type in reaction {orig_idx}."
                    )
                plog = True
            else:
                plog = False
            falloff = reaction.rate.type == "falloff"
            is_troe = reaction.rate.sub_type == "Troe"
            is_sri = reaction.rate.sub_type == "Sri"
            is_lindemann = reaction.rate.sub_type == "Lindemann"
            aeuc = cu.activation_energy_units()
            if not third_body and not falloff and not plog:
                # Case 3 !PD, !TB, !PLOG
                ctuc = cu.prefactor_units(cc.ureg("kmol/m**3"), 1 - dim)
                pef = (reaction.rate.pre_exponential_factor * ctuc).to_base_units()
                beta = reaction.rate.temperature_exponent
                ae = (
                    reaction.rate.activation_energy * cc.ureg.joule / cc.ureg.kmol
                ).to(aeuc)
            elif not falloff and not plog:
                # Case 2 !PD, TB, !PLOG
                ctuc = cu.prefactor_units(cc.ureg("kmol/m**3"), -dim)
                pef = (reaction.rate.pre_exponential_factor * ctuc).to_base_units()
                beta = reaction.rate.temperature_exponent
                ae = (
                    reaction.rate.activation_energy * cc.ureg.joule / cc.ureg.kmol
                ).to(aeuc)
            elif not third_body and not falloff and plog:
                # Case 4 PLOG
                ctuc = cu.prefactor_units(cc.ureg("kmol/m**3"), 1 - dim)
                plog_pef, plog_beta, plog_ae = cu.evaluate_plog(
                    reaction.rate.rates, mechanism.P
                )
                pef = (plog_pef * ctuc).to_base_units()
                beta = plog_beta
                ae = (plog_ae * cc.ureg.joule / cc.ureg.kmol).to(aeuc)
            elif third_body and not falloff and plog:
                # Case 5 third body PLOG
                ctuc = cu.prefactor_units(cc.ureg("kmol/m**3"), -dim)
                plog_pef, plog_beta, plog_ae = cu.evaluate_plog(
                    reaction.rate.rates, mechanism.P
                )
                pef = (plog_pef * ctuc).to_base_units()
                beta = plog_beta
                ae = (plog_ae * cc.ureg.joule / cc.ureg.kmol).to(aeuc)
            else:
                # Case 1 PD, TB
                ctuc = cu.prefactor_units(cc.ureg("kmol/m**3"), 1 - dim)
                pef = (
                    reaction.rate.high_rate.pre_exponential_factor * ctuc
                ).to_base_units()
                beta = reaction.rate.high_rate.temperature_exponent
                ae = (
                    reaction.rate.high_rate.activation_energy
                    * cc.ureg.joule
                    / cc.ureg.kmol
                ).to(aeuc)

                low_pef = (
                    reaction.rate.low_rate.pre_exponential_factor * ctuc
                ).to_base_units()
                low_beta = reaction.rate.low_rate.temperature_exponent
                low_ae = (
                    reaction.rate.low_rate.activation_energy
                    * cc.ureg.joule
                    / cc.ureg.kmol
                ).to(aeuc)
                if is_troe:
                    troe = reaction.rate.falloff_coeffs
                    ntroe = len(troe)
                elif is_sri:
                    sri = reaction.rate.falloff_coeffs
                    nsri = len(sri)
                elif is_lindemann:
                    pass
                else:
                    raise ValueError(
                        f"Unrecognized reaction rate type {reaction.rate.type},"
                        f" {reaction.rate.sub_type} for reaction: {reaction.equation}"
                    )

            cw.writer(
                fstream,
                cw.comment(f"reaction {orig_idx}:  {reaction.equation}"),
            )
            cw.writer(fstream, f"k_f = {pef.m:.15g}")
            if (beta == 0) and (ae == 0):
                cw.writer(fstream, "           ;")
            else:
                if ae == 0:
                    cw.writer(fstream, f"           * exp(({beta:.15g}) * logT);")
                elif beta == 0:
                    cw.writer(
                        fstream,
                        "           *"
                        f" exp(-({(1.0 / cc.Rc / cc.ureg.kelvin * ae).m:.15g})"
                        " * invT);",
                    )
                else:
                    cw.writer(
                        fstream,
                        f"           * exp(({beta:.15g}) * logT -"
                        f" ({(1.0 / cc.Rc / cc.ureg.kelvin * ae).m:.15g})"
                        " * invT);",
                    )

            alpha = None
            if not third_body and not falloff:
                cw.writer(fstream, f"qf[{idx}] *= k_f;")
            elif (
                not falloff
                and len(reaction.third_body.efficiencies) == 1
                and isclose(reaction.third_body.default_efficiency, 0.0)
            ):
                cw.writer(fstream, f"qf[{idx}] *= k_f;")
            elif not falloff:
                alpha = enhancement_d_with_qss(mechanism, species_info, reaction)
                cw.writer(fstream, f"Corr  = {alpha};")
                cw.writer(fstream, f"qf[{idx}] *= Corr * k_f;")
            else:
                alpha = enhancement_d_with_qss(mechanism, species_info, reaction)
                cw.writer(fstream, f"Corr  = {alpha};")
                cw.writer(
                    fstream,
                    "redP = Corr / k_f *"
                    f" {10 ** (-dim * 6) * low_pef.m * 10 ** 3 ** dim:.15g} ",
                )
                if (low_beta == 0) and (low_ae.m == 0):
                    cw.writer(
                        fstream,
                        "           ;",
                    )
                elif low_ae.m == 0:
                    cw.writer(
                        fstream,
                        f"           * exp({low_beta:.15g}  * logT);",
                    )
                elif low_beta == 0:
                    cw.writer(
                        fstream,
                        "           * exp(-"
                        f" ({(1.0 / cc.Rc / cc.ureg.kelvin * low_ae).m:.15g})"
                        " *invT);",
                    )
                else:
                    cw.writer(
                        fstream,
                        f"           * exp({low_beta:.15g}  * logT -"
                        f" ({(1.0 / cc.Rc / cc.ureg.kelvin * low_ae).m:.15g})"
                        " *invT);",
                    )
                if is_troe:
                    cw.writer(fstream, "F = redP / (1.0 + redP);")
                    cw.writer(fstream, "logPred = log10(redP);")
                    cw.writer(fstream, "logFcent = log10(")
                    if abs(troe[1]) > 1.0e-100:
                        if 1.0 - troe[0] != 0:
                            cw.writer(
                                fstream,
                                f"    ({1.0 - troe[0]:.15g})*exp(-T *"
                                f" {1 / troe[1]:.15g})",
                            )
                    else:
                        cw.writer(fstream, "     0.0 ")
                    if abs(troe[2]) > 1.0e-100:
                        if troe[0] == 1:
                            cw.writer(
                                fstream,
                                f"    + exp(-T * {1 / troe[2]:.15g})",
                            )
                        else:
                            cw.writer(
                                fstream,
                                f"    + {troe[0]:.15g} * exp(-T * {1 / troe[2]:.15g})",
                            )
                    else:
                        cw.writer(fstream, "    + 0.0 ")
                    if ntroe == 4:
                        if troe[3] < 0:
                            cw.writer(fstream, f"    + exp({-troe[3]:.15g} * invT));")
                        else:
                            cw.writer(fstream, f"    + exp(-{troe[3]:.15g} * invT));")
                    else:
                        cw.writer(fstream, "    + 0.0);")
                    cw.writer(fstream, "troe_c = -0.4 - 0.67 * logFcent;")
                    cw.writer(fstream, "troe_n = 0.75 - 1.27 * logFcent;")
                    cw.writer(
                        fstream,
                        "troe = (troe_c + logPred) / (troe_n -"
                        " 0.14*(troe_c + logPred));",
                    )
                    cw.writer(
                        fstream,
                        "F_troe = exp(M_LN10 * logFcent / (1.0 + troe*troe));",
                    )
                    cw.writer(fstream, "Corr = F * F_troe;")
                    cw.writer(fstream, f"qf[{idx}] *= Corr * k_f;")
                elif is_sri:
                    cw.writer(fstream, "F = redP / (1.0 + redP);")
                    cw.writer(fstream, "logPred = log10(redP);")
                    cw.writer(fstream, "X = 1.0 / (1.0 + logPred*logPred);")
                    if sri[1] < 0:
                        cw.writer(
                            fstream,
                            f"F_sri = exp(X * log({sri[0]:.15g} *"
                            f" exp({-sri[1]:.15g}*invT)",
                        )
                    else:
                        cw.writer(
                            fstream,
                            f"F_sri = exp(X * log({sri[0]:.15g} *"
                            f" exp(-{sri[1]:.15g}*invT)",
                        )
                    if sri[2] > 1.0e-100:
                        cw.writer(fstream, f"   +  exp(logT/{sri[2]:.15g}) ")
                    else:
                        cw.writer(fstream, "   +  0. ")
                    cw.writer(
                        fstream,
                        f"   *  ({nsri} > 3 ?"
                        f" {sri[3]:.15g}*exp({sri[4]:.15g}*logT : 1.0);",
                    )
                    cw.writer(fstream, "Corr = F * F_sri;")
                    cw.writer(fstream, f"qf[{idx}] *= Corr * k_f;")
                elif nlindemann > 0:
                    cw.writer(fstream, "Corr = redP / (1. + redP);")
                    cw.writer(fstream, f"qf[{idx}] *= Corr * k_f;")

            if kc_conv_inv:
                if alpha is None:
                    cw.writer(
                        fstream,
                        f"qr[{idx}] *= k_f * {kc_e} * ({kc_conv_inv});",
                    )
                else:
                    cw.writer(
                        fstream,
                        f"qr[{idx}] *= Corr * k_f * {kc_e} *" f" ({kc_conv_inv});",
                    )
            else:
                if alpha is None:
                    cw.writer(
                        fstream,
                        f"qr[{idx}] *= k_f * {kc_e};",
                    )
                else:
                    cw.writer(
                        fstream,
                        f"qr[{idx}] *= Corr * k_f * {kc_e};",
                    )

        cw.writer(fstream)

    cw.writer(fstream)
    cw.writer(fstream, "}")

    cw.writer(fstream)

    # main function
    if n_reactions > 0:
        cw.writer(
            fstream,
            "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void "
            " productionRate(amrex::Real * wdot, const amrex::Real * sc,"
            " const amrex::Real T)",
        )
    else:
        cw.writer(
            fstream,
            "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void "
            " productionRate(amrex::Real * wdot, const amrex::Real *"
            " /*sc*/, const amrex::Real /*T*/)",
        )
    cw.writer(fstream, "{")

    if n_reactions == 0:
        cw.writer(fstream)
    else:
        cw.writer(fstream, "const amrex::Real invT = 1.0 / T;")
        cw.writer(fstream, "const amrex::Real logT = log(T);")
        cw.writer(fstream)
        cw.writer(
            fstream,
            cw.comment("reference concentration: P_atm / (RT) in inverse mol/m^3"),
        )
        cw.writer(
            fstream,
            f"const amrex::Real refC = {cc.Patm_pa:g} /"
            f" {cc.R.to(cc.ureg.joule / (cc.ureg.mole / cc.ureg.kelvin)).m:g} *"
            " invT;",
        )
        cw.writer(fstream, "const amrex::Real refCinv = 1 / refC;")

        if nsri > 0:
            cw.writer(fstream, "amrex::Real X, F_sri;")

    cw.writer(fstream)
    cw.writer(fstream, f"for (int i = 0; i < {n_species}; ++i) {{")
    cw.writer(fstream, "wdot[i] = 0.0;")
    cw.writer(fstream, "}")
    cw.writer(fstream)

    # initialize the symbolic wdot array
    for i in range(n_species):
        syms.wdot_smp[i] = 0.0

    # initialize the symbolic jacobian array
    for i in range(n_species * (n_species + 1)):
        syms.jac_smp[i] = 0.0

    if n_reactions > 0:
        # nclassd = n_reactions - nspecial
        # nCorr   = n3body + ntroe + nsri + nlindemann

        # Mixt concentration for PD & TB
        cw.writer(fstream, cw.comment("compute the mixture concentration"))
        cw.writer(fstream, "amrex::Real mixture = 0.0;")
        cw.writer(fstream, f"for (int i = 0; i < {n_species}; ++i) {{")
        cw.writer(fstream, "mixture += sc[i];")
        cw.writer(fstream, "}")
        cw.writer(fstream)

        # Kc stuff
        cw.writer(fstream, cw.comment("compute the Gibbs free energy"))
        cw.writer(fstream, f"amrex::Real g_RT[{species_info.n_species}];")
        cw.writer(fstream, "gibbs(g_RT, T);")
        if species_info.n_qssa_species > 0:
            cw.writer(
                fstream,
                f"amrex::Real g_RT_qss[{species_info.n_qssa_species}];",
            )
            cw.writer(fstream, "gibbs_qss(g_RT_qss, T);")
        cu.write_gibbs_exp(fstream, mechanism, species_info)
        cw.writer(fstream)

        if species_info.n_qssa_species > 0:
            cw.writer(
                fstream,
                f"amrex::Real sc_qss[{max(1, species_info.n_qssa_species)}];",
            )
            cw.writer(
                fstream,
                f"amrex::Real kf_qss[{reaction_info.n_qssa_reactions}],"
                f" qf_qss[{reaction_info.n_qssa_reactions}],"
                f" qr_qss[{reaction_info.n_qssa_reactions}];",
            )
            cw.writer(fstream, cw.comment("Fill sc_qss here"))
            cw.writer(fstream, "comp_k_f_qss(T, invT, logT, kf_qss);")
            # cw.writer(fstream,"comp_Kc_qss(invT, g_RT, g_RT_qss, Kc_qss);")
            cw.writer(
                fstream,
                "comp_qss_coeff(kf_qss, qf_qss, qr_qss, sc, T, g_RT, g_RT_qss);",
            )
            cw.writer(fstream, "comp_sc_qss(sc_qss, qf_qss, qr_qss);")
            cw.writer(fstream)

        # Loop like you're going through them in the mech.Linp order
        for orig_idx, _ in reaction_info.idxmap.items():
            reaction = mechanism.reaction(orig_idx)
            cw.writer(fstream, "{")
            if bool(reaction.orders):
                forward_sc, forward_sc_smp = cu.qss_sorted_phase_space(
                    mechanism, species_info, reaction, reaction.orders, syms
                )
            else:
                forward_sc, forward_sc_smp = cu.qss_sorted_phase_space(
                    mechanism, species_info, reaction, reaction.reactants, syms
                )
            if reaction.reversible:
                reverse_sc, reverse_sc_smp = cu.qss_sorted_phase_space(
                    mechanism, species_info, reaction, reaction.products, syms
                )
            else:
                reverse_sc = "0.0"
                reverse_sc_smp = 0.0

            kc_exp_arg, kc_exp_arg_smp = cu.sorted_kc_exp_arg(
                mechanism, species_info, reaction, syms
            )

            kc_e = cu.kc_exp(mechanism, species_info, reaction, -1)
            kc_conv_inv, kc_conv_inv_smp = cu.fkc_conv_inv(
                mechanism, species_info, reaction, syms
            )

            if bool(reaction.orders):
                dim = cu.phase_space_units(reaction.orders)
            else:
                dim = cu.phase_space_units(reaction.reactants)

            third_body = reaction.third_body is not None
            if reaction.reaction_type == "pressure-dependent-Arrhenius":
                plog = True
            elif reaction.reaction_type == "three-body-pressure-dependent-Arrhenius":
                if not third_body:
                    raise ValueError(
                        f"Inconsistent reaction type in reaction {orig_idx}."
                    )
                plog = True
            else:
                plog = False
            falloff = reaction.rate.type == "falloff"
            is_troe = reaction.rate.sub_type == "Troe"
            is_sri = reaction.rate.sub_type == "Sri"
            is_lindemann = reaction.rate.sub_type == "Lindemann"
            aeuc = cu.activation_energy_units()
            if not third_body and not falloff and not plog:
                # Case 3 !PD, !TB
                ctuc = cu.prefactor_units(cc.ureg("kmol/m**3"), 1 - dim)
                pef = (reaction.rate.pre_exponential_factor * ctuc).to_base_units()
                beta = reaction.rate.temperature_exponent
                ae = (
                    reaction.rate.activation_energy * cc.ureg.joule / cc.ureg.kmol
                ).to(aeuc)
            elif not falloff and not plog:
                # Case 2 !PD, TB, !PLOG
                ctuc = cu.prefactor_units(cc.ureg("kmol/m**3"), -dim)
                pef = (reaction.rate.pre_exponential_factor * ctuc).to_base_units()
                beta = reaction.rate.temperature_exponent
                ae = (
                    reaction.rate.activation_energy * cc.ureg.joule / cc.ureg.kmol
                ).to(aeuc)
            elif not third_body and not falloff and plog:
                # Case 4 PLOG
                ctuc = cu.prefactor_units(cc.ureg("kmol/m**3"), 1 - dim)
                plog_pef, plog_beta, plog_ae = cu.evaluate_plog(
                    reaction.rate.rates, mechanism.P
                )
                pef = (plog_pef * ctuc).to_base_units()
                beta = plog_beta
                ae = (plog_ae * cc.ureg.joule / cc.ureg.kmol).to(aeuc)
            elif third_body and not falloff and plog:
                # Case 5 third body PLOG
                ctuc = cu.prefactor_units(cc.ureg("kmol/m**3"), -dim)
                plog_pef, plog_beta, plog_ae = cu.evaluate_plog(
                    reaction.rate.rates, mechanism.P
                )
                pef = (plog_pef * ctuc).to_base_units()
                beta = plog_beta
                ae = (plog_ae * cc.ureg.joule / cc.ureg.kmol).to(aeuc)
            else:
                # Case 1 PD, TB, !PLOG
                ctuc = cu.prefactor_units(cc.ureg("kmol/m**3"), 1 - dim)
                pef = (
                    reaction.rate.high_rate.pre_exponential_factor * ctuc
                ).to_base_units()
                beta = reaction.rate.high_rate.temperature_exponent
                ae = (
                    reaction.rate.high_rate.activation_energy
                    * cc.ureg.joule
                    / cc.ureg.kmol
                ).to(aeuc)

                low_pef = (
                    reaction.rate.low_rate.pre_exponential_factor * ctuc
                ).to_base_units()
                low_beta = reaction.rate.low_rate.temperature_exponent
                low_ae = (
                    reaction.rate.low_rate.activation_energy
                    * cc.ureg.joule
                    / cc.ureg.kmol
                ).to(aeuc)
                if is_troe:
                    troe = reaction.rate.falloff_coeffs
                    ntroe = len(troe)
                elif is_sri:
                    sri = reaction.rate.falloff_coeffs
                    nsri = len(sri)
                elif is_lindemann:
                    pass
                else:
                    raise ValueError(
                        f"Unrecognized reaction rate type {reaction.rate.type},"
                        f" {reaction.rate.sub_type} for reaction: {reaction.equation}"
                    )

                low_beta = syms.convert_number_to_int(low_beta)

            beta = syms.convert_number_to_int(beta)

            cw.writer(
                fstream,
                cw.comment(f"reaction {orig_idx}:  {reaction.equation}"),
            )
            cw.writer(fstream, f"const amrex::Real k_f = {pef.m:.15g}")

            k_f_smp = pef.m

            if (beta == 0) and (ae == 0):
                cw.writer(fstream, "           ;")
            else:
                if ae == 0:
                    cw.writer(fstream, f"           * exp(({beta:.15g}) * logT);")
                    k_f_smp *= sme.exp(beta * syms.logT_smp)
                elif beta == 0:
                    cw.writer(
                        fstream,
                        "           *"
                        f" exp(-({(1.0 / cc.Rc / cc.ureg.kelvin * ae).m:.15g})"
                        " * invT);",
                    )
                    coeff = (((1.0 / cc.Rc / cc.ureg.kelvin)) * ae).magnitude
                    k_f_smp *= sme.exp(-coeff * syms.invT_smp)
                else:
                    cw.writer(
                        fstream,
                        f"           * exp(({beta:.15g}) * logT -"
                        f" ({(1.0 / cc.Rc / cc.ureg.kelvin * ae).m:.15g})"
                        " * invT);",
                    )
                    coeff = ((1.0 / cc.Rc / cc.ureg.kelvin)) * ae
                    k_f_smp *= sme.exp(beta * syms.logT_smp - coeff * syms.invT_smp)

            alpha = None
            if not third_body and not falloff:
                cw.writer(
                    fstream,
                    f"const amrex::Real qf = k_f * ({forward_sc});",
                )
                qf_smp = k_f_smp * forward_sc_smp
            elif (
                not falloff
                and len(reaction.third_body.efficiencies) == 1
                and isclose(reaction.third_body.default_efficiency, 0.0)
            ):
                cw.writer(
                    fstream,
                    f"const amrex::Real qf = k_f * ({forward_sc});",
                )
                qf_smp = k_f_smp * forward_sc_smp
            elif not falloff:
                alpha, alpha_smp = enhancement_d_with_qss(
                    mechanism, species_info, reaction, syms
                )
                cw.writer(fstream, f"const amrex::Real Corr = {alpha};")
                corr_smp = alpha_smp
                cw.writer(
                    fstream,
                    f"const amrex::Real qf = Corr * k_f * ({forward_sc});",
                )
                qf_smp = corr_smp * k_f_smp * forward_sc_smp
            else:
                alpha, alpha_smp = enhancement_d_with_qss(
                    mechanism, species_info, reaction, syms
                )
                cw.writer(fstream, f"amrex::Real Corr = {alpha};")
                corr_smp = alpha_smp
                cw.writer(
                    fstream,
                    "const amrex::Real redP = Corr / k_f *"
                    f" {10 ** (-dim * 6) * low_pef.m * 10 ** 3 ** dim:.15g} ",
                )
                redp_smp = (
                    corr_smp / k_f_smp * (10 ** (-dim * 6) * low_pef.m * 10 ** (3**dim))
                )
                if (low_beta == 0) and (low_ae.m == 0):
                    cw.writer(
                        fstream,
                        "           ;",
                    )
                elif low_ae.m == 0:
                    cw.writer(
                        fstream,
                        f"           * exp({low_beta:.15g} * logT);",
                    )
                elif low_beta == 0:
                    cw.writer(
                        fstream,
                        "           * exp(-"
                        f" {(1.0 / cc.Rc / cc.ureg.kelvin * low_ae).m:.15g} *"
                        " invT);",
                    )
                else:
                    cw.writer(
                        fstream,
                        f"           * exp({low_beta:.15g} * logT -"
                        f" {(1.0 / cc.Rc / cc.ureg.kelvin * low_ae).m:.15g} *"
                        " invT);",
                    )
                coeff = (1.0 / cc.Rc / cc.ureg.kelvin * low_ae).magnitude
                redp_smp *= sme.exp(low_beta * syms.logT_smp - coeff * syms.invT_smp)
                if is_troe:
                    cw.writer(fstream, "const amrex::Real F = redP / (1.0 + redP);")
                    f_smp = redp_smp / (1.0 + redp_smp)
                    cw.writer(fstream, "const amrex::Real logPred = log10(redP);")
                    logpred_smp = sme.log(redp_smp, 10)
                    cw.writer(fstream, "const amrex::Real logFcent = log10(")
                    int_smp = 0.0
                    if abs(troe[1]) > 1.0e-100:
                        if 1.0 - troe[0] != 0:
                            cw.writer(
                                fstream,
                                f"    {1.0 - troe[0]:.15g} * exp(-T *"
                                f" {1 / troe[1]:.15g})",
                            )
                            first_factor = syms.convert_number_to_int(1.0 - troe[0])
                            second_factor = syms.convert_number_to_int(-1 / troe[1])
                            int_smp += first_factor * sme.exp(
                                syms.T_smp * second_factor
                            )
                    else:
                        cw.writer(fstream, "     0.0 ")
                    if abs(troe[2]) > 1.0e-100:
                        if troe[0] != 0:
                            cw.writer(
                                fstream,
                                f"    + {troe[0]:.15g} * exp(-T * {1 / troe[2]:.15g})",
                            )
                            first_factor = syms.convert_number_to_int(troe[0])
                            second_factor = syms.convert_number_to_int(-1 / troe[2])
                            int_smp += first_factor * sme.exp(
                                syms.T_smp * second_factor
                            )
                    else:
                        cw.writer(fstream, "    + 0.0 ")
                    if ntroe == 4:
                        if troe[3] < 0:
                            cw.writer(fstream, f"    + exp({-troe[3]:.15g} * invT));")
                            first_factor = syms.convert_number_to_int(-troe[3])
                            int_smp += sme.exp(first_factor * syms.invT_smp)
                        else:
                            cw.writer(fstream, f"    + exp(-{troe[3]:.15g} * invT));")
                            first_factor = syms.convert_number_to_int(-troe[3])
                            int_smp += sme.exp(first_factor * syms.invT_smp)
                    else:
                        cw.writer(fstream, "    + 0.0);")
                    logfcent_smp = sme.log(int_smp, 10)
                    cw.writer(
                        fstream,
                        "const amrex::Real troe_c = -0.4 - 0.67 * logFcent;",
                    )
                    troe_c_smp = -0.4 - 0.67 * logfcent_smp
                    cw.writer(
                        fstream,
                        "const amrex::Real troe_n = 0.75 - 1.27 * logFcent;",
                    )
                    troe_n_smp = 0.75 - 1.27 * logfcent_smp
                    cw.writer(
                        fstream,
                        "const amrex::Real troe = (troe_c + logPred) /"
                        " (troe_n - 0.14 * (troe_c + logPred));",
                    )
                    troe_smp = (troe_c_smp + logpred_smp) / (
                        troe_n_smp - 0.14 * (troe_c_smp + logpred_smp)
                    )
                    cw.writer(
                        fstream,
                        "const amrex::Real F_troe = exp(M_LN10 *logFcent /"
                        " (1.0 + troe * troe));",
                    )
                    f_troe_smp = pow(10, logfcent_smp / (1.0 + troe_smp * troe_smp))
                    cw.writer(fstream, "Corr = F * F_troe;")
                    corr_smp = f_smp * f_troe_smp
                    cw.writer(
                        fstream,
                        f"const amrex::Real qf = Corr * k_f * ({forward_sc});",
                    )
                    qf_smp = corr_smp * k_f_smp * forward_sc_smp
                elif is_sri:
                    cw.writer(fstream, "const amrex::Real F = redP / (1.0 + redP);")
                    f_smp = redp_smp / (1.0 + redp_smp)
                    cw.writer(fstream, "const amrex::Real logPred = log10(redP);")
                    logpred_smp = sme.log(redp_smp, 10)
                    cw.writer(fstream, "X = 1.0 / (1.0 + logPred*logPred);")
                    # x_smp = 1.0 / (1.0 + logpred_smp * logpred_smp)
                    if sri[1] < 0:
                        cw.writer(
                            fstream,
                            f"F_sri = exp(X * log({sri[0]:.15g} *"
                            f" exp({-sri[1]:.15g} * invT)",
                        )
                        if syms is not None:
                            raise NotImplementedError("Not done for now")
                    else:
                        cw.writer(
                            fstream,
                            f"F_sri = exp(X * log({sri[0]:.15g} *"
                            f" exp(-{sri[1]:.15g} * invT)",
                        )
                        if syms is not None:
                            raise NotImplementedError("Not done for now")
                    if sri[2] > 1.0e-100:
                        cw.writer(fstream, f"   +  exp(logT / {sri[2]:.15g}) ")
                        if syms is not None:
                            raise NotImplementedError("Not done for now")
                    else:
                        cw.writer(fstream, "   +  0. ")
                        if syms is not None:
                            raise NotImplementedError("Not done for now")
                    cw.writer(
                        fstream,
                        f"   *  ({nsri} > 3 ? {sri[3]:.15g} *"
                        f" exp({sri[4]:.15g} * logT) : 1.0);",
                    )
                    if syms is not None:
                        raise NotImplementedError("Not done for now")
                    cw.writer(fstream, "Corr = F * F_sri;")
                    cw.writer(
                        fstream,
                        f"const amrex::Real qf = Corr * k_f * ({forward_sc});",
                    )
                elif nlindemann > 0:
                    cw.writer(fstream, "Corr = redP / (1.0 + redP);")
                    corr_smp = redp_smp / (1.0 + redp_smp)
                    cw.writer(
                        fstream,
                        f"const amrex::Real qf = Corr * k_f * ({forward_sc});",
                    )
                    qf_smp = corr_smp * k_f_smp * forward_sc_smp
            if kc_conv_inv:
                if alpha is None:
                    if reverse_sc_smp == 0:
                        cw.writer(
                            fstream,
                            "const amrex::Real qr = 0.0;",
                        )
                    else:
                        cw.writer(
                            fstream,
                            "const amrex::Real qr = k_f *"
                            f" {kc_e} * ({kc_conv_inv}) *"
                            f" ({reverse_sc});",
                        )
                    qr_smp = (
                        k_f_smp
                        * sme.exp(-(kc_exp_arg_smp))
                        * (kc_conv_inv_smp)
                        * (reverse_sc_smp)
                    )
                else:
                    cw.writer(
                        fstream,
                        "const amrex::Real qr = Corr * k_f *"
                        f" {kc_e} * ({kc_conv_inv}) *"
                        f" ({reverse_sc});",
                    )
                    qr_smp = (
                        corr_smp
                        * k_f_smp
                        * sme.exp(-kc_exp_arg_smp)
                        * (kc_conv_inv_smp)
                        * (reverse_sc_smp)
                    )
            else:
                if alpha is None:
                    if reverse_sc_smp == 0.0:
                        cw.writer(
                            fstream,
                            "const amrex::Real qr = 0.0;",
                        )
                    else:
                        cw.writer(
                            fstream,
                            "const amrex::Real qr = k_f *" f" {kc_e} * ({reverse_sc});",
                        )
                    qr_smp = k_f_smp * sme.exp(-(kc_exp_arg_smp)) * reverse_sc_smp
                else:
                    cw.writer(
                        fstream,
                        "const amrex::Real qr = Corr * k_f *"
                        f" {kc_e} * ({reverse_sc});",
                    )
                    qr_smp = (
                        corr_smp * k_f_smp * sme.exp(-(kc_exp_arg_smp)) * reverse_sc_smp
                    )

            remove_forward = cu.is_remove_forward(reaction_info, orig_idx)

            if remove_forward:
                cw.writer(fstream, cw.comment("Remove forward reaction"))
                cw.writer(fstream, cw.comment("const amrex::Real qdot = qf - qr;"))
                cw.writer(fstream, "const amrex::Real qdot = - qr;")
                qdot_smp = -qr_smp
            else:
                cw.writer(fstream, "const amrex::Real qdot = qf - qr;")
                qdot_smp = qf_smp - qr_smp

            reaction = mechanism.reaction(orig_idx)
            lst_reactants = [(k, v) for k, v in reaction.reactants.items()]
            lst_products = [(k, v) for k, v in reaction.products.items()]
            all_agents = list(set(lst_reactants + lst_products))
            agents = []
            # remove QSS species from agents
            for symbol, coefficient in all_agents:
                if symbol not in species_info.qssa_species_list:
                    agents.append((symbol, coefficient))
            dict_species = {v: i for i, v in enumerate(species_info.all_species_list)}
            agents = sorted(
                agents,
                key=lambda v, dict_species=dict_species: dict_species[v[0]],
            )
            # Check for duplicates
            assert len(agents) == len(
                set(agents)
            ), f"Reaction {reaction} contains duplicate agents"
            # note that a species might appear as both reactant and product
            # a species might also appear twice or more on on each side
            # agents is a set that contains unique (symbol, coefficient)
            for a in agents:
                symbol, coefficient = a
                for b in reaction.reactants:
                    if b == a[0] and reaction.reactants[b] == a[1]:
                        if coefficient == 1.0:
                            cw.writer(
                                fstream,
                                f"wdot[{species_info.ordered_idx_map[symbol]}]"
                                " -= qdot;",
                            )
                            syms.wdot_smp[
                                species_info.ordered_idx_map[symbol]
                            ] -= qdot_smp
                        else:
                            cw.writer(
                                fstream,
                                f"wdot[{species_info.ordered_idx_map[symbol]}]"
                                f" -= {coefficient:f} * qdot;",
                            )
                            coefficient = syms.convert_number_to_int(coefficient)
                            syms.wdot_smp[species_info.ordered_idx_map[symbol]] -= (
                                coefficient * qdot_smp
                            )
                for b in reaction.products:
                    if b == a[0] and reaction.products[b] == a[1]:
                        if coefficient == 1.0:
                            cw.writer(
                                fstream,
                                f"wdot[{species_info.ordered_idx_map[symbol]}]"
                                " += qdot;",
                            )
                            syms.wdot_smp[
                                species_info.ordered_idx_map[symbol]
                            ] += qdot_smp
                        else:
                            cw.writer(
                                fstream,
                                f"wdot[{species_info.ordered_idx_map[symbol]}]"
                                f" += {coefficient:f} * qdot;",
                            )
                            coefficient = syms.convert_number_to_int(coefficient)
                            syms.wdot_smp[species_info.ordered_idx_map[symbol]] += (
                                coefficient * qdot_smp
                            )

            cw.writer(fstream, "}")
            cw.writer(fstream)

    cw.writer(fstream)

    cw.writer(fstream, "}")

    cw.writer(fstream)


def production_rate_light(fstream, mechanism, species_info, reaction_info):
    """Write low memory production rate."""
    n_species = species_info.n_species
    n_reactions = mechanism.n_reactions

    assert len(reaction_info.index) == 7

    itroe = reaction_info.index[0:2]
    isri = reaction_info.index[1:3]
    ilindemann = reaction_info.index[2:4]
    # i3body = reaction_info.index[3:5]
    # isimple = reaction_info.index[4:6]
    # ispecial = reaction_info.index[5:7]

    ntroe = itroe[1] - itroe[0]
    nsri = isri[1] - isri[0]
    nlindemann = ilindemann[1] - ilindemann[0]
    # n3body = i3body[1] - i3body[0]
    # nsimple = isimple[1] - isimple[0]
    # nspecial = ispecial[1] - ispecial[0]

    # qdot
    cw.writer(fstream)
    cw.writer(
        fstream,
        cw.comment(
            "GPU version of productionRate: no more use of thermo namespace vectors "
        ),
    )
    # main function
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void "
        "productionRate_light("
        "amrex::Real * wdot,"
        "amrex::Real * sc,"
        "amrex::Real * g_RT,"
        "amrex::Real * g_RT_qss,"
        "amrex::Real * sc_qss,"
        "amrex::Real * kf_qss,"
        "amrex::Real * qf_qss,"
        "amrex::Real * qr_qss,"
        "const amrex::Real T,"
        "const amrex::Real invT,"
        "const amrex::Real logT)",
    )
    cw.writer(fstream, "{")

    if n_reactions == 0:
        cw.writer(fstream)
    else:
        cw.writer(
            fstream,
            cw.comment("reference concentration: P_atm / (RT) in inverse mol/m^3"),
        )
        cw.writer(
            fstream,
            f"const amrex::Real refC = {cc.Patm_pa:g} /"
            f" {cc.R.to(cc.ureg.joule / (cc.ureg.mole / cc.ureg.kelvin)).m:g} *"
            " invT;",
        )
        cw.writer(fstream, "const amrex::Real refCinv = 1 / refC;")

        if nsri > 0:
            cw.writer(fstream, "amrex::Real X, F_sri;")

    cw.writer(fstream)
    cw.writer(fstream, f"for (int i = 0; i < {n_species}; ++i) {{")
    cw.writer(fstream, "wdot[i] = 0.0;")
    cw.writer(fstream, "}")
    cw.writer(fstream)

    if n_reactions > 0:
        # nclassd = n_reactions - nspecial
        # nCorr   = n3body + ntroe + nsri + nlindemann

        # Mixt concentration for PD & TB
        cw.writer(fstream, cw.comment("compute the mixture concentration"))
        cw.writer(fstream, "amrex::Real mixture = 0.0;")
        cw.writer(fstream, f"for (int i = 0; i < {n_species}; ++i) {{")
        cw.writer(fstream, "mixture += sc[i];")
        cw.writer(fstream, "}")
        cw.writer(fstream)

        # Kc stuff
        cw.writer(fstream, cw.comment("compute the Gibbs free energy"))
        cw.writer(fstream, "gibbs(g_RT, T);")
        if species_info.n_qssa_species > 0:
            cw.writer(fstream, "gibbs_qss(g_RT_qss, T);")
        cu.write_gibbs_exp(fstream, mechanism, species_info)
        cw.writer(fstream)

        if species_info.n_qssa_species > 0:
            cw.writer(fstream, cw.comment("Fill sc_qss here"))
            cw.writer(fstream, "comp_k_f_qss(T, invT, logT, kf_qss);")
            # cw.writer(fstream,"comp_Kc_qss(invT, g_RT, g_RT_qss, Kc_qss);")
            cw.writer(
                fstream,
                "comp_qss_coeff(kf_qss, qf_qss, qr_qss, sc, T, g_RT, g_RT_qss);",
            )
            cw.writer(fstream, "comp_sc_qss(sc_qss, qf_qss, qr_qss);")
            cw.writer(fstream)

        # Loop like you're going through them in the mech.Linp order
        for orig_idx, _ in reaction_info.idxmap.items():
            reaction = mechanism.reaction(orig_idx)
            cw.writer(fstream, "{")
            if bool(reaction.orders):
                forward_sc = cu.qss_sorted_phase_space(
                    mechanism, species_info, reaction, reaction.orders
                )
            else:
                forward_sc = cu.qss_sorted_phase_space(
                    mechanism, species_info, reaction, reaction.reactants
                )
            if reaction.reversible:
                reverse_sc = cu.qss_sorted_phase_space(
                    mechanism, species_info, reaction, reaction.products
                )
            else:
                reverse_sc = "0.0"

            kc_e = cu.kc_exp(mechanism, species_info, reaction, -1)
            kc_conv_inv = cu.fkc_conv_inv(mechanism, species_info, reaction)

            if bool(reaction.orders):
                dim = cu.phase_space_units(reaction.orders)
            else:
                dim = cu.phase_space_units(reaction.reactants)
            third_body = reaction.third_body is not None
            falloff = reaction.rate.type == "falloff"
            is_troe = reaction.rate.sub_type == "Troe"
            is_sri = reaction.rate.sub_type == "Sri"
            is_lindemann = reaction.rate.sub_type == "Lindemann"
            aeuc = cu.activation_energy_units()
            if not third_body and not falloff:
                # Case 3 !PD, !TB
                ctuc = cu.prefactor_units(cc.ureg("kmol/m**3"), 1 - dim)
                pef = (reaction.rate.pre_exponential_factor * ctuc).to_base_units()
                beta = reaction.rate.temperature_exponent
                ae = (
                    reaction.rate.activation_energy * cc.ureg.joule / cc.ureg.kmol
                ).to(aeuc)
            elif not falloff:
                # Case 2 !PD, TB
                ctuc = cu.prefactor_units(cc.ureg("kmol/m**3"), -dim)
                pef = (reaction.rate.pre_exponential_factor * ctuc).to_base_units()
                beta = reaction.rate.temperature_exponent
                ae = (
                    reaction.rate.activation_energy * cc.ureg.joule / cc.ureg.kmol
                ).to(aeuc)
            else:
                # Case 1 PD, TB
                ctuc = cu.prefactor_units(cc.ureg("kmol/m**3"), 1 - dim)
                pef = (
                    reaction.rate.high_rate.pre_exponential_factor * ctuc
                ).to_base_units()
                beta = reaction.rate.high_rate.temperature_exponent
                ae = (
                    reaction.rate.high_rate.activation_energy
                    * cc.ureg.joule
                    / cc.ureg.kmol
                ).to(aeuc)

                low_pef = (
                    reaction.rate.low_rate.pre_exponential_factor * ctuc
                ).to_base_units()
                low_beta = reaction.rate.low_rate.temperature_exponent
                low_ae = (
                    reaction.rate.low_rate.activation_energy
                    * cc.ureg.joule
                    / cc.ureg.kmol
                ).to(aeuc)
                if is_troe:
                    troe = reaction.rate.falloff_coeffs
                    ntroe = len(troe)
                    is_troe = True
                elif is_sri:
                    sri = reaction.rate.falloff_coeffs
                    nsri = len(sri)
                elif is_lindemann:
                    pass
                else:
                    raise ValueError(
                        f"Unrecognized reaction rate type {reaction.rate.type},"
                        f" {reaction.rate.sub_type} for reaction: {reaction.equation}"
                    )

            cw.writer(
                fstream,
                cw.comment(f"reaction {orig_idx}:  {reaction.equation}"),
            )
            cw.writer(fstream, f"const amrex::Real k_f = {pef.m:.15g}")

            if (beta == 0) and (ae == 0):
                cw.writer(fstream, "           ;")
            else:
                if ae == 0:
                    cw.writer(fstream, f"           * exp(({beta:.15g}) * logT);")
                elif beta == 0:
                    cw.writer(
                        fstream,
                        "           *"
                        f" exp(-({(1.0 / cc.Rc / cc.ureg.kelvin * ae).m:.15g})"
                        " * invT);",
                    )
                else:
                    cw.writer(
                        fstream,
                        f"           * exp(({beta:.15g}) * logT -"
                        f" ({(1.0 / cc.Rc / cc.ureg.kelvin * ae).m:.15g}) *"
                        " invT);",
                    )

            alpha = None
            if not third_body and not falloff:
                cw.writer(
                    fstream,
                    f"const amrex::Real qf = k_f * ({forward_sc});",
                )
            elif (
                not falloff
                and len(reaction.third_body.efficiencies) == 1
                and isclose(reaction.third_body.default_efficiency, 0.0)
            ):
                cw.writer(
                    fstream,
                    f"const amrex::Real qf = k_f * ({forward_sc});",
                )
            elif not falloff:
                alpha = enhancement_d_with_qss(mechanism, species_info, reaction)
                cw.writer(fstream, f"const amrex::Real Corr = {alpha};")
                cw.writer(
                    fstream,
                    f"const amrex::Real qf = Corr * k_f * ({forward_sc});",
                )
            else:
                alpha = enhancement_d_with_qss(mechanism, species_info, reaction)
                cw.writer(fstream, f"amrex::Real Corr = {alpha};")
                cw.writer(
                    fstream,
                    "const amrex::Real redP = Corr / k_f *"
                    f" {10 ** (-dim * 6) * low_pef.m * 10 ** 3 ** dim:.15g} ",
                )
                cw.writer(
                    fstream,
                    f"           * exp({low_beta:.15g} * logT -"
                    f" {(1.0 / cc.Rc / cc.ureg.kelvin * low_ae).m:.15g} *"
                    " invT);",
                )
                if is_troe:
                    cw.writer(fstream, "const amrex::Real F = redP / (1.0 + redP);")
                    cw.writer(fstream, "const amrex::Real logPred = log10(redP);")
                    cw.writer(fstream, "const amrex::Real logFcent = log10(")
                    if abs(troe[1]) > 1.0e-100:
                        if 1.0 - troe[0] != 0:
                            cw.writer(
                                fstream,
                                f"    {1.0 - troe[0]:.15g} * exp(-T *"
                                f" {1 / troe[1]:.15g})",
                            )
                    else:
                        cw.writer(fstream, "     0.0 ")
                    if abs(troe[2]) > 1.0e-100:
                        if troe[0] != 0:
                            cw.writer(
                                fstream,
                                f"    + {troe[0]:.15g} * exp(-T * {1 / troe[2]:.15g})",
                            )
                    else:
                        cw.writer(fstream, "    + 0.0 ")
                    if ntroe == 4:
                        if troe[3] < 0:
                            cw.writer(fstream, f"    + exp({-troe[3]:.15g} * invT));")
                        else:
                            cw.writer(fstream, f"    + exp(-{troe[3]:.15g} * invT));")
                    else:
                        cw.writer(fstream, "    + 0.0);")
                    cw.writer(
                        fstream,
                        "const amrex::Real troe_c = -0.4 - 0.67 * logFcent;",
                    )
                    cw.writer(
                        fstream,
                        "const amrex::Real troe_n = 0.75 - 1.27 * logFcent;",
                    )
                    cw.writer(
                        fstream,
                        "const amrex::Real troe = (troe_c + logPred) /"
                        " (troe_n - 0.14 * (troe_c + logPred));",
                    )
                    cw.writer(
                        fstream,
                        "const amrex::Real F_troe = exp(M_LN10 *logFcent /"
                        " (1.0 + troe * troe));",
                    )
                    cw.writer(fstream, "Corr = F * F_troe;")
                    cw.writer(
                        fstream,
                        f"const amrex::Real qf = Corr * k_f * ({forward_sc});",
                    )
                elif is_sri:
                    cw.writer(fstream, "const amrex::Real F = redP / (1.0 + redP);")
                    cw.writer(fstream, "const amrex::Real logPred = log10(redP);")
                    cw.writer(fstream, "X = 1.0 / (1.0 + logPred*logPred);")
                    if sri[1] < 0:
                        cw.writer(
                            fstream,
                            f"F_sri = exp(X * log({sri[0]:.15g} *"
                            f" exp({-sri[1]:.15g} * invT)",
                        )
                    else:
                        cw.writer(
                            fstream,
                            f"F_sri = exp(X * log({sri[0]:.15g} *"
                            f" exp(-{sri[1]:.15g} * invT)",
                        )
                    if sri[2] > 1.0e-100:
                        cw.writer(fstream, f"   +  exp(logT / {sri[2]:.15g}) ")
                    else:
                        cw.writer(fstream, "   +  0. ")
                    cw.writer(
                        fstream,
                        f"   *  ({nsri} > 3 ? {sri[3]:.15g} *"
                        f" exp({sri[4]:.15g} * logT) : 1.0);",
                    )
                    cw.writer(fstream, "Corr = F * F_sri;")
                    cw.writer(
                        fstream,
                        f"const amrex::Real qf = Corr * k_f * ({forward_sc});",
                    )
                elif nlindemann > 0:
                    cw.writer(fstream, "Corr = redP / (1.0 + redP);")
                    cw.writer(
                        fstream,
                        f"const amrex::Real qf = Corr * k_f * ({forward_sc});",
                    )
            if kc_conv_inv:
                if alpha is None:
                    if reverse_sc == "0.0":
                        cw.writer(
                            fstream,
                            "const amrex::Real qr = 0.0;",
                        )
                    else:
                        cw.writer(
                            fstream,
                            "const amrex::Real qr = k_f *"
                            f" {kc_e} * ({kc_conv_inv}) *"
                            f" ({reverse_sc});",
                        )
                else:
                    cw.writer(
                        fstream,
                        "const amrex::Real qr = Corr * k_f *"
                        f" {kc_e} * ({kc_conv_inv}) *"
                        f" ({reverse_sc});",
                    )
            else:
                if alpha is None:
                    if reverse_sc == "0.0":
                        cw.writer(
                            fstream,
                            "const amrex::Real qr = 0.0;",
                        )
                    else:
                        cw.writer(
                            fstream,
                            "const amrex::Real qr = k_f *" f" {kc_e} * ({reverse_sc});",
                        )
                else:
                    cw.writer(
                        fstream,
                        "const amrex::Real qr = Corr * k_f *"
                        f" {kc_e} * ({reverse_sc});",
                    )
            remove_forward = cu.is_remove_forward(reaction_info, orig_idx)

            if remove_forward:
                cw.writer(fstream, cw.comment("Remove forward reaction"))
                cw.writer(fstream, cw.comment("const amrex::Real qdot = qf - qr;"))
                cw.writer(fstream, "const amrex::Real qdot = - qr;")
            else:
                cw.writer(fstream, "const amrex::Real qdot = qf - qr;")

            reaction = mechanism.reaction(orig_idx)
            lst_reactants = [(k, v) for k, v in reaction.reactants.items()]
            lst_products = [(k, v) for k, v in reaction.products.items()]
            all_agents = list(set(lst_reactants + lst_products))
            agents = []
            # remove QSS species from agents
            for symbol, coefficient in all_agents:
                if symbol not in species_info.qssa_species_list:
                    agents.append((symbol, coefficient))
            dict_species = {v: i for i, v in enumerate(species_info.all_species_list)}
            agents = sorted(
                agents,
                key=lambda v, dict_species=dict_species: dict_species[v[0]],
            )
            # Check for duplicates
            assert len(agents) == len(
                set(agents)
            ), f"Reaction {reaction} contains duplicate agents"
            # note that a species might appear as both reactant and product
            # a species might also appear twice or more on on each side
            # agents is a set that contains unique (symbol, coefficient)
            for a in agents:
                symbol, coefficient = a
                for b in reaction.reactants:
                    if b == a[0] and reaction.reactants[b] == a[1]:
                        if coefficient == 1.0:
                            cw.writer(
                                fstream,
                                f"wdot[{species_info.ordered_idx_map[symbol]}]"
                                " -= qdot;",
                            )
                        else:
                            cw.writer(
                                fstream,
                                f"wdot[{species_info.ordered_idx_map[symbol]}]"
                                f" -= {coefficient:f} * qdot;",
                            )
                for b in reaction.products:
                    if b == a[0] and reaction.products[b] == a[1]:
                        if coefficient == 1.0:
                            cw.writer(
                                fstream,
                                f"wdot[{species_info.ordered_idx_map[symbol]}]"
                                " += qdot;",
                            )
                        else:
                            cw.writer(
                                fstream,
                                f"wdot[{species_info.ordered_idx_map[symbol]}]"
                                f" += {coefficient:f} * qdot;",
                            )
            cw.writer(fstream, "}")
            cw.writer(fstream)

    cw.writer(fstream)

    cw.writer(fstream, "return;")
    cw.writer(fstream, "}")

    cw.writer(fstream)


def progress_rate_fr(fstream, species_info, reaction_info):
    """Write progress rates."""
    assert len(reaction_info.index) == 7

    cw.writer(fstream)
    cw.writer(fstream, cw.comment("compute the progress rate for each reaction"))
    cw.writer(fstream, cw.comment("USES progressRate : todo switch to GPU"))
    if reaction_info.n_reactions > 0:
        cw.writer(
            fstream,
            "void progressRateFR"
            + "(amrex::Real *  q_f, amrex::Real *  q_r, amrex::Real *  sc,"
            " amrex::Real T)",
        )
    else:
        cw.writer(
            fstream,
            "void progressRateFR"
            + "(amrex::Real *  /*q_f*/, amrex::Real *  /*q_r*/, amrex::Real *  /*sc*/,"
            " amrex::Real /*T*/)",
        )

    cw.writer(fstream, "{")

    if reaction_info.n_reactions > 0:
        cw.writer(fstream, "const amrex::Real invT = 1.0 / T;")
        cw.writer(fstream, "const amrex::Real logT = log(T);")

        cw.writer(fstream, cw.comment("compute the Gibbs free energy"))
        cw.writer(fstream, "amrex::Real g_RT[NUM_SPECIES];")
        cw.writer(fstream, "gibbs(g_RT, T);")
        if species_info.n_qssa_species > 0:
            cw.writer(
                fstream,
                "amrex::Real g_RT_qss[NUM_QSS_GAS_SPECIES];",
            )
            cw.writer(fstream, "gibbs_qss(g_RT_qss, T);")

        cw.writer(fstream)
        cw.writer(
            fstream,
            f"amrex::Real sc_qss[{max(1, species_info.n_qssa_species)}];",
        )
        if species_info.n_qssa_species > 0:
            cw.writer(fstream, cw.comment("Fill sc_qss here"))
            cw.writer(
                fstream,
                f"amrex::Real kf_qss[{reaction_info.n_qssa_reactions}],"
                f" qf_qss[{reaction_info.n_qssa_reactions}],"
                f" qr_qss[{reaction_info.n_qssa_reactions}];",
            )
            cw.writer(fstream, "comp_k_f_qss(T, invT, logT, kf_qss);")
            cw.writer(
                fstream,
                "comp_qss_coeff(kf_qss, qf_qss, qr_qss, sc, T, g_RT, g_RT_qss);",
            )
            cw.writer(fstream, "comp_sc_qss(sc_qss, qf_qss, qr_qss);")

        cw.writer(fstream, "comp_qfqr(q_f, q_r, sc, sc_qss, T, invT, logT);")
        cw.writer(fstream)

    cw.writer(fstream, "}")


def enhancement_d_with_qss(mechanism, species_info, reaction, syms=None):
    """Write enhancement with QSS."""
    record_symbolic_operations = True
    if syms is None:
        record_symbolic_operations = False
    third_body = reaction.third_body is not None
    falloff = reaction.rate.type == "falloff"
    if not third_body and not falloff:
        raise ValueError("enhancement_d called for a reaction without a third body")

    if not reaction.third_body:
        raise NotImplementedError("FIXME EFFICIENCIES")
        species, coefficient = third_body
        if species == "<mixture>":
            if record_symbolic_operations:
                return "mixture", syms.mixture_smp
            else:
                return "mixture"
        if record_symbolic_operations:
            return (
                f"sc[{species_info.ordered_idx_map[species]}]",
                syms.sc_smp[species_info.ordered_idx_map[species]],
            )
        else:
            return f"sc[{species_info.ordered_idx_map[species]}]"

    efficiencies = reaction.third_body.efficiencies
    alpha = ["mixture"]
    if record_symbolic_operations:
        alpha_smp = [syms.mixture_smp]
    dict_species = {v: i for i, v in enumerate(species_info.all_species_list)}
    sorted_efficiencies = sorted(efficiencies.keys(), key=lambda v: dict_species[v])
    for symbol in sorted_efficiencies:
        efficiency = efficiencies[symbol]
        if symbol not in species_info.qssa_species_list:
            factor = f"({efficiency - 1:.15g})"
            if record_symbolic_operations:
                factor_smp = efficiency - 1
            if (efficiency - 1) != 0:
                conc = f"sc[{species_info.ordered_idx_map[symbol]}]"
                if record_symbolic_operations:
                    conc_smp = syms.sc_smp[species_info.ordered_idx_map[symbol]]
                if (efficiency - 1) == 1:
                    alpha.append(f"{conc}")
                    if record_symbolic_operations:
                        alpha_smp.append(conc_smp)
                else:
                    alpha.append(f"{factor}*{conc}")
                    if record_symbolic_operations:
                        alpha_smp.append(factor_smp * conc_smp)
        else:
            factor = f"({efficiency - 1:.15g})"
            if record_symbolic_operations:
                factor_smp = efficiency - 1
            if (efficiency - 1) != 0:
                idx = species_info.ordered_idx_map[symbol] - species_info.n_species
                conc = f"sc_qss[{idx}]"
                if record_symbolic_operations:
                    conc_smp = syms.sc_qss_smp[
                        species_info.ordered_idx_map[symbol] - species_info.n_species
                    ]
                if (efficiency - 1) == 1:
                    alpha.append(f"{conc}")
                    if record_symbolic_operations:
                        alpha_smp.append(conc_smp)
                else:
                    alpha.append(f"{factor}*{conc}")
                    if record_symbolic_operations:
                        alpha_smp.append(factor_smp * conc_smp)

    if record_symbolic_operations:
        enhancement_smp = 0.0
        for alpha_val in alpha_smp:
            enhancement_smp += alpha_val

    if not record_symbolic_operations:
        return " + ".join(alpha).replace("+ -", "- ")
    else:
        return " + ".join(alpha).replace("+ -", "- "), enhancement_smp


def surface_chemistry(fstream, interface, species_info):
    """Write the surface reaction routines of a heterogeneous mechanism.

    Surface reactions are numbered in mechanism order, as in SKINU.
    """
    site_density = interface.site_density * 1000.0  # kmol/m^2 -> mol/m^2
    # surface standard concentrations site_density / sigma_k [mol/m^2]
    ref_conc = {sp.name: site_density / sp.size for sp in interface.species()}
    reactions = interface.reactions()
    surface_progress_rate(fstream, reactions, species_info, site_density, ref_conc)
    surface_production_rate(fstream, reactions, species_info)
    surface_kfkr(fstream, species_info, ref_conc)


def _int_power(base, exponent):
    """Return the C++ product base^exponent for a positive integer exponent."""
    return " * ".join([base] * int(exponent))


def _concentration_product(species_info, orders):
    """Return the C++ product of concentrations sc[k]^order (sc in SI mol units).

    As in Cantera, a negative order on an absent species gives a zero rate.
    """
    terms = []
    for name, order in orders.items():
        sc = f"sc[{species_info.ordered_idx_map[name]}]"
        if order > 0.0 and float(order).is_integer():
            terms.append(_int_power(sc, order))
        elif order > 0.0:
            terms.append(f"pow({sc}, {order:.15g})")
        elif order < 0.0:
            terms.append(f"({sc} > 0.0 ? pow({sc}, {order:.15g}) : 0.0)")
    return " * ".join(terms) if terms else "1.0"


def _coverage_terms(species_info, reaction, ref_conc):
    """Return the C++ exponent terms of a surface reaction's coverage dependencies.

    A dependency on species k multiplies the rate constant by
    exp(a theta_k ln10 + m ln(max(theta_k, 1e-20)) - E theta_k / RT), as in
    Cantera, with the coverage theta_k = sc[k] / ref_conc_k. The clamp keeps
    negative m finite on a bare surface.
    """
    terms = []
    for name, dep in reaction.rate.coverage_dependencies.items():
        idx = species_info.ordered_idx_map[name]
        theta = f"(sc[{idx}] * {1.0 / ref_conc[name]:.15g})"
        if dep["a"] != 0.0:
            terms.append(f"{dep['a'] * log(10.0):+.15g} * {theta}")
        if dep["m"] != 0.0:
            terms.append(f"{dep['m']:+.15g} * log({theta} > 1e-20 ? {theta} : 1e-20)")
        if dep["E"] != 0.0:
            # E is in J/kmol
            terms.append(f"{-dep['E'] / ct.gas_constant:+.15g} * {theta} * invT")
    return terms


def _arrhenius(pef, beta, ea_over_r, extra_terms=()):
    """Return a C++ modified Arrhenius expression with extra exponent terms."""
    arg = []
    if beta != 0.0:
        arg.append(f"{beta:+.15g} * logT")
    if ea_over_r != 0.0:
        arg.append(f"{-ea_over_r:+.15g} * invT")
    arg += extra_terms
    if not arg:
        return f"{pef:.15g}"
    return f"{pef:.15g} * exp({' '.join(arg).lstrip('+')})"


def surface_progress_rate(fstream, reactions, species_info, site_density, ref_conc):
    """Write comp_qfqr_surface: forward/reverse rates of the surface reactions.

    sc (NUM_ALL_SPECIES) holds gas concentrations [mol/m^3] at
    0 .. NUM_GAS_SPECIES-1, followed by surface concentrations [mol/m^2]
    (theta_k * site_density / sigma_k).
    qf/qr (NUM_SURFACE_REACTIONS, in mechanism order as in SKINU) are
    returned in [mol/m^2/s]. Rates follow Cantera's interface and sticking
    Arrhenius forms, including coverage dependencies, reaction orders and the
    Motz-Wise correction. Reverse rates use k_r = k_f / Kc, with standard
    concentrations P_atm/RT (gas) and site_density/sigma_k (surface).
    """
    r_kmol = ct.gas_constant  # J/kmol/K
    # gas standard concentration P_atm / (RT) = ref_c_gas / T in mol/m^3
    ref_c_gas = cc.Patm_pa / (r_kmol / 1000.0)
    reversible_species = {
        name for r in reactions if r.reversible for name in cu.kc_net_coefficients(r)
    }
    need_surface_g = bool(reversible_species & ref_conc.keys())
    need_gas_g = bool(reversible_species - ref_conc.keys())
    any_sticking = any(r.reaction_type == "sticking-Arrhenius" for r in reactions)

    cw.writer(fstream)
    cw.writer(fstream, cw.comment("Forward and reverse rates of the surface reactions"))
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void comp_qfqr_surface(amrex::Real *"
        " qf, amrex::Real * qr, const amrex::Real * sc, const amrex::Real T, const"
        " amrex::Real invT, const amrex::Real logT)",
    )
    cw.writer(fstream, "{")
    cw.writer(fstream, "amrex::ignore_unused(T, invT, logT);")
    if need_gas_g:
        cw.writer(fstream, "amrex::Real g_RT[NUM_GAS_SPECIES];")
        cw.writer(fstream, "gibbs(g_RT, T);")
    if need_surface_g:
        cw.writer(fstream, "amrex::Real g_RT_s[NUM_SURFACE_SPECIES];")
        cw.writer(fstream, "gibbs_surface(g_RT_s, T);")
    if any_sticking:
        cw.writer(fstream, "const amrex::Real sqrtT = sqrt(T);")
    cw.writer(fstream)

    for i, reaction in enumerate(reactions):
        rate = reaction.rate
        orders = dict(reaction.reactants)
        orders.update(reaction.orders)
        ea_over_r = rate.activation_energy / r_kmol
        coverage = _coverage_terms(species_info, reaction, ref_conc)

        cw.writer(fstream, cw.comment(f"surface reaction {i}: {reaction.equation}"))
        cw.writer(fstream, "{")
        if reaction.reaction_type == "sticking-Arrhenius":
            # gamma is dimensionless:
            # k_f = gamma / site_density^m * sqrt(RT / (2 pi W))
            gamma = _arrhenius(
                rate.pre_exponential_factor,
                rate.temperature_exponent,
                ea_over_r,
                coverage,
            )
            cw.writer(fstream, f"amrex::Real gamma = {gamma};")
            if rate.motz_wise_correction:
                cw.writer(fstream, "gamma /= (1.0 - 0.5 * gamma);")
            weight = rate.sticking_weight / 1000.0  # kg/mol
            coeff = sqrt(r_kmol / 1000.0 / (2.0 * pi * weight)) / (
                site_density**rate.sticking_order
            )
            k_f = f"gamma * {coeff:.15g} * sqrtT"
        else:
            # interface-Arrhenius (others are rejected by check_surface_reactions)
            # kmol-based SI prefactor -> mol-based: 1000^(1 - total order)
            pef = rate.pre_exponential_factor * 1000.0 ** (
                1.0 - cu.phase_space_units(orders)
            )
            k_f = _arrhenius(pef, rate.temperature_exponent, ea_over_r, coverage)
        cw.writer(fstream, f"const amrex::Real k_f = {k_f};")

        cw.writer(
            fstream,
            f"qf[{i}] = k_f * {_concentration_product(species_info, orders)};",
        )

        if reaction.reversible:
            # 1/Kc = exp(sum_k nu_k g_k/RT) * prod_k (C0_k)^-nu_k, with
            # nu_k = products - reactants
            dg = []
            conv = 1.0
            dn_gas = 0.0
            for name, minus_n in cu.kc_net_coefficients(reaction).items():
                n = -minus_n
                if name in ref_conc:
                    k = species_info.ordered_idx_map[name] - species_info.n_gas_species
                    dg.append(f"{n:+.15g} * g_RT_s[{k}]")
                    conv *= ref_conc[name] ** (-n)
                else:
                    k = species_info.ordered_idx_map[name]
                    dg.append(f"{n:+.15g} * g_RT[{k}]")
                    dn_gas += n
            kc_inv = f"exp({' '.join(dg)}) * {conv:.15g}"
            if dn_gas != 0.0:
                if float(dn_gas).is_integer():
                    # (ref_c_gas / T)^-dn_gas
                    base = (
                        f"({1.0 / ref_c_gas:.15g} * T)"
                        if dn_gas > 0
                        else f"({ref_c_gas:.15g} * invT)"
                    )
                    kc_inv += f" * {_int_power(base, abs(dn_gas))}"
                else:
                    kc_inv += f" * pow({ref_c_gas:.15g} * invT, {-dn_gas:.15g})"
            products = _concentration_product(species_info, reaction.products)
            cw.writer(fstream, f"qr[{i}] = k_f * {kc_inv} * {products};")
        else:
            cw.writer(fstream, f"qr[{i}] = 0.0;")
        cw.writer(fstream, "}")
    cw.writer(fstream, "}")


def surface_production_rate(fstream, reactions, species_info):
    """Write surfaceProductionRate: per-area production rates of all species.

    sdot [mol/m^2/s] covers all NUM_ALL_SPECIES (gas, then surface) and only
    the surface reactions; sc is laid out as for comp_qfqr_surface.
    """
    cw.writer(fstream)
    cw.writer(
        fstream,
        cw.comment("Production rates of all species from surface reactions"),
    )
    cw.writer(fstream, cw.comment("sdot in mol/m^2/s; sc as in comp_qfqr_surface"))
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void surfaceProductionRate("
        "amrex::Real * sdot, const amrex::Real * sc, const amrex::Real T)",
    )
    cw.writer(fstream, "{")
    cw.writer(fstream, "const amrex::Real invT = 1.0 / T;")
    cw.writer(fstream, "const amrex::Real logT = log(T);")
    cw.writer(
        fstream,
        "amrex::Real qf[NUM_SURFACE_REACTIONS], qr[NUM_SURFACE_REACTIONS];",
    )
    cw.writer(fstream, "comp_qfqr_surface(qf, qr, sc, T, invT, logT);")
    cw.writer(fstream, "for (int i = 0; i < NUM_ALL_SPECIES; ++i) {")
    cw.writer(fstream, "sdot[i] = 0.0;")
    cw.writer(fstream, "}")
    for i, reaction in enumerate(reactions):
        cw.writer(fstream, "{")
        cw.writer(fstream, f"const amrex::Real qdot = qf[{i}] - qr[{i}];")
        for name, coef in reaction.reactants.items():
            k = species_info.ordered_idx_map[name]
            cw.writer(fstream, f"sdot[{k}] -= {coef:.15g} * qdot;")
        for name, coef in reaction.products.items():
            k = species_info.ordered_idx_map[name]
            cw.writer(fstream, f"sdot[{k}] += {coef:.15g} * qdot;")
        cw.writer(fstream, "}")
    cw.writer(fstream, "}")


def surface_kfkr(fstream, species_info, ref_conc):
    """Write SKKFKR: surface progress rates from P, T, mole fractions, coverages.

    x (NUM_ALL_SPECIES) holds gas mole fractions followed by surface
    coverages; q_f/q_r (NUM_SURFACE_REACTIONS) are returned in chemkin units
    [mol/cm^2/s]. P is in dyn/cm^2.
    """
    cw.writer(fstream)
    cw.writer(
        fstream, cw.comment("Returns the progress rates of each surface reaction")
    )
    cw.writer(
        fstream,
        cw.comment("Given P, T, gas mole fractions and surface coverages"),
    )
    cw.writer(
        fstream,
        "AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void SKKFKR(const amrex::Real P,"
        " const amrex::Real T, const amrex::Real x[], amrex::Real q_f[],"
        " amrex::Real q_r[])",
    )
    cw.writer(fstream, "{")
    cw.writer(fstream, "amrex::Real sc[NUM_ALL_SPECIES];")
    cw.writer(
        fstream,
        "const amrex::Real PORT = 1e6 *"
        f" P/({(cc.R * cc.ureg.kelvin * cc.ureg.mole / cc.ureg.erg).m:1.14e} *"
        " T); " + cw.comment("1e6 * P/RT so sc goes to SI units"),
    )
    cw.writer(fstream, "for (int id = 0; id < NUM_GAS_SPECIES; ++id) {")
    cw.writer(fstream, "sc[id] = x[id] * PORT;")
    cw.writer(fstream, "}")
    cw.writer(
        fstream,
        cw.comment("surface concentrations: theta_k * site_density / sigma_k"),
    )
    for name in species_info.surface_species_list:
        idx = species_info.ordered_idx_map[name]
        cw.writer(
            fstream,
            f"sc[{idx}] = x[{idx}] * {ref_conc[name]:.15g};" + cw.comment(name),
        )
    cw.writer(
        fstream,
        "comp_qfqr_surface(q_f, q_r, sc, T, 1.0 / T, log(T));",
    )
    cw.writer(fstream, cw.comment("convert mol/m^2/s to chemkin units, mol/cm^2/s"))
    cw.writer(fstream, "for (int id = 0; id < NUM_SURFACE_REACTIONS; ++id) {")
    cw.writer(fstream, "q_f[id] *= 1.0e-4;")
    cw.writer(fstream, "q_r[id] *= 1.0e-4;")
    cw.writer(fstream, "}")
    cw.writer(fstream, "}")
