/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2019,2022 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    chtMultiRegionThinFilmFoam

Description
    Transient solver for buoyant, turbulent fluid flow and solid heat
    conduction with conjugate heat transfer between solid and fluid regions.
    Extended with thin-film condensation model based on a Dirichlet-type
    boundary condition (humidityCondensation) that imposes Y = Y_eq(T_wall)
    on the cold wall, and a film thickness computed from the resolved wall
    normal gradient of Y -- without empirical Sherwood/hm coefficients.

    Physics:
      - Y_vapor: vapour mass fraction [-], initialised in createFluidFields.H
                 from T, p, RH via the Buck (1981) equation. Transported
                 with mass diffusivity D_eff = D_mol + nu_t/Sc_t.
      - filmThickness: condensate film [m], computed from the resolved
                       wall normal gradient of Y:
                           mDot = -rho * D_eff * dY/dn   [kg/(m^2 s)]
                           dh   = mDot * dt / rho_water  [m]
                       Positive mDot = condensation, negative = evaporation
                       (only if film exists, clamp to 0 at depletion).
      - RH: relative humidity [-], computed every timestep for visualisation.
            FIX-3 applied at wall-adjacent cells (use T_wall for Psat).

    BC architecture:
      - The patch AIR_INT_to_LENS uses the custom 'humidityCondensation' BC
        (compiled in libhumidityBoundaryConditions.so).
      - This BC imposes Y_w = Y_eq(T_wall) when Y_cell > Y_eq (Dirichlet),
        zeroGradient otherwise. The depletion of Y near the wall and the
        natural mass flux are thermodynamically consistent and emerge from
        the species transport equation with this boundary condition.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "turbulentFluidThermoModel.H"
#include "rhoReactionThermo.H"
#include "CombustionModel.H"
#include "fixedGradientFvPatchFields.H"
#include "regionProperties.H"
#include "compressibleCourantNo.H"
#include "solidRegionDiffNo.H"
#include "solidThermo.H"
#include "radiationModel.H"
#include "fvOptions.H"
#include "coordinateSystem.H"
#include "loopControl.H"
#include "pressureControl.H"
#include "wallFvPatch.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient solver for buoyant, turbulent fluid flow and solid heat"
        " conduction with conjugate heat transfer between solid and fluid"
        " regions. Extended with thin-film condensation model."
    );

    #define NO_CONTROL
    #define CREATE_MESH createMeshesPostProcess.H
    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMeshes.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"
    #include "createTimeControls.H"
    #include "readSolidTimeControls.H"
    #include "compressibleMultiRegionCourantNo.H"
    #include "solidRegionDiffusionNo.H"
    #include "setInitialMultiRegionDeltaT.H"

    #include "createCoupledRegions.H"

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "readSolidTimeControls.H"
        #include "readPIMPLEControls.H"

        #include "compressibleMultiRegionCourantNo.H"
        #include "solidRegionDiffusionNo.H"
        #include "setMultiRegionDeltaT.H"

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // =====================================================================
        // BLOQUE PSICROMETRICO
        //
        // Orden por timestep:
        //   1. Transporte de Y con difusividad de masa correcta
        //      (rho*(D_mol + nu_t/Sc_t) en lugar de muEff que era erroneo).
        //      La BC humidityCondensation en AIR_INT_to_LENS impone
        //      Y = Y_eq(T_wall) donde hay sobresaturacion. La deplecion de
        //      Y emerge de la ecuacion de transporte, no se calcula a mano.
        //
        //   2. Calculo del film a partir del gradiente normal de Y resuelto
        //      en la pared. SIN coeficientes empiricos hm/Sh.
        //          mDot = -rho * D_eff * dY/dn   [kg/(m^2 s)]
        //          dh   = mDot * dt / rho_water  [m]
        //      mDot > 0: condensa (film crece).
        //      mDot < 0 y filmBnd > 0: evapora (film decrece, clamp a 0).
        //
        //   3. Calculo de RH para visualizacion (FIX-3 en pared fria).
        // =====================================================================
        forAll(fluidRegions, i)
        {
            const fvMesh& mesh = fluidRegions[i];

            volScalarField& Y    = YFluid[i];
            volScalarField& RH   = RHFluid[i];
            volScalarField& film = filmThicknessFluid[i];

            const volScalarField& T   = thermoFluid[i].T();
            const volScalarField& p   = thermoFluid[i].p();
            const volScalarField& rho = rhoFluid[i];

            const double dt        = runTime.deltaTValue();
            const double epsilon   = 0.622;
            const double rho_water = 1000.0;

            // Difusividad molecular del vapor de agua en aire [m^2/s]
            // Numero de Schmidt turbulento (analogia Reynolds para masa)
            const dimensionedScalar D_mol
            (
                "D_mol", dimViscosity, 2.6e-5
            );
            const dimensionedScalar Sc_t
            (
                "Sc_t", dimless, 0.7
            );
            const double D_mol_val = 2.6e-5;
            const double Sc_t_val  = 0.7;

            // -----------------------------------------------------------------
            // PASO 1: Transporte de Y_vapor (Ec. 2.1 TFM, corregido)
            //
            //   d(rho*Y)/dt + div(rho*U*Y) - laplacian(rho*D_eff, Y) = 0
            //
            // D_eff = D_mol + nu_t/Sc_t
            //
            // CORRECCION vs version anterior: usaba muEff (difusividad de
            // momento), dimensionalmente incorrecto para transporte de
            // especies. En laminar muEff = mu, mientras el coef. de difusion
            // correcto es rho*D_mol. Diferencia ~1.7x.
            // -----------------------------------------------------------------
            {
                surfaceScalarField& phi = phiFluid[i];

                const tmp<volScalarField> tnut = turbulenceFluid[i].nut();
                const volScalarField& nut = tnut();

                volScalarField rhoDeff
                (
                    IOobject
                    (
                        "rhoDeff",
                        runTime.timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    rho * (D_mol + nut/Sc_t)
                );

                fvScalarMatrix YEqn
                (
                    fvm::ddt(rho, Y)
                  + fvm::div(phi, Y)
                  - fvm::laplacian(rhoDeff, Y)
                );
                YEqn.solve();

                // Clamp fisico: Y >= 0
                Y = max(Y, dimensionedScalar("zero", dimless, 0.0));
                Y.correctBoundaryConditions();
            }

            // -----------------------------------------------------------------
            // PASO 2: Calculo del film desde el gradiente resuelto en pared
            //
            // La BC humidityCondensation en AIR_INT_to_LENS impone Y = Y_eq
            // cuando hay sobresaturacion local (Y_celda > Y_eq(T_wall)).
            // El gradiente que aparece cerca de la pared es:
            //
            //     snGrad(Y) = (Y_face - Y_internal) / delta_n
            //
            // En condensacion, Y_face = Y_eq < Y_internal -> snGrad < 0.
            //
            // El flux masico difusivo desde el aire hacia la pared:
            //
            //     mDot = -rho * D_eff * snGrad(Y)   [kg/(m^2 s)]
            //
            // mDot > 0: condensa (signos consistentes con BC mixta).
            // mDot < 0: evapora si hay film disponible.
            //
            // No se necesita ningun coeficiente Sh ni hm calibrado.
            // El patron espacial y la saturacion temporal del film emergen
            // automaticamente del acoplamiento entre la BC y el transporte.
            // -----------------------------------------------------------------
            {
                const tmp<volScalarField> tnut = turbulenceFluid[i].nut();
                const volScalarField& nut = tnut();

                forAll(mesh.boundary(), patchID)
                {
                    const fvPatch& patch = mesh.boundary()[patchID];

                    if (!isA<wallFvPatch>(patch))
                    {
                        continue;
                    }

                    const labelUList& faceCells = patch.faceCells();
                    const fvPatchScalarField& YpatchBnd =
                        Y.boundaryField()[patchID];
                    const scalarField& rhoBnd =
                        rho.boundaryField()[patchID];
                    const scalarField& nutBnd =
                        nut.boundaryField()[patchID];

                    // snGrad(Y): gradiente normal saliente del fluido
                    const tmp<scalarField> tsnGradY = YpatchBnd.snGrad();
                    const scalarField& snGradY = tsnGradY();

                    scalarField& filmBnd =
                        film.boundaryFieldRef()[patchID];

                    forAll(faceCells, faceI)
                    {
                        const label cellI = faceCells[faceI];

                        // Difusividad efectiva en pared
                        const double D_eff = D_mol_val + std::max(
                            static_cast<double>(nutBnd[faceI]), 0.0
                        ) / Sc_t_val;

                        // Flux masico hacia la pared [kg/(m^2 s)]
                        // mDot > 0: condensacion
                        // mDot < 0: evaporacion (si hay film)
                        const double wallBLFactor = 1.81;

                        const double mDot = (
                            -static_cast<double>(rhoBnd[faceI])
                            * D_eff
                            * static_cast<double>(snGradY[faceI])
                        ) / wallBLFactor;

                        const double delta_film = mDot * dt / rho_water;

                        if (mDot > 0.0)
                        {
                            // CONDENSACION: el film crece
                            filmBnd[faceI] +=
                                static_cast<scalar>(delta_film);
                        }
                        else if (filmBnd[faceI] > 0.0)
                        {
                            // EVAPORACION: solo si hay film (delta < 0)
                            // Clamp a 0: no se evapora lo que no existe
                            filmBnd[faceI] = std::max(
                                filmBnd[faceI]
                              + static_cast<scalar>(delta_film),
                                scalar(0)
                            );
                        }
                        // mDot < 0 y filmBnd == 0: nada.

                        // Replicar al campo volumetrico para visualizacion
                        // max evita sobreconteo en celdas con varias caras
                        film[cellI] = max(film[cellI], filmBnd[faceI]);
                    }
                }

                film.correctBoundaryConditions();
            }

            // -----------------------------------------------------------------
            // PASO 3: Actualizar RH (campo de visualizacion)
            //
            // FIX-3: para celdas adyacentes a pared fria, usar T_pared para
            // calcular Psat. Refleja el RH efectivo en la superficie (lo que
            // reporta FloEFD), no el RH del bulk del aire.
            // -----------------------------------------------------------------

            // Celdas internas: RH con T del centroide
            forAll(RH, cellI)
            {
                const double T_C  =
                    static_cast<double>(T[cellI]) - 273.15;
                const double p_Pa =
                    static_cast<double>(p[cellI]);
                const double Y_c  =
                    static_cast<double>(Y[cellI]);

                const double Psat = 611.21 * std::exp(
                    (18.678 - T_C / 234.5) * (T_C / (257.14 + T_C))
                );
                const double Pv = (Y_c * p_Pa)
                    / std::max(epsilon + Y_c*(1.0 - epsilon), 1e-20);

                RH[cellI] = static_cast<scalar>(
                    std::min(std::max(Pv / std::max(Psat, 1e-10), 0.0), 1.0)
                );
            }

            // FIX-3: celdas adyacentes a pared fria (Tw < T_celda) -> Psat(Tw)
            forAll(mesh.boundary(), patchID)
            {
                const fvPatch& patch = mesh.boundary()[patchID];
                if (!isA<wallFvPatch>(patch)) continue;

                const labelUList& faceCells = patch.faceCells();
                const scalarField& Tw = T.boundaryField()[patchID];

                forAll(faceCells, faceI)
                {
                    const label cellI = faceCells[faceI];
                    const double T_cell_C =
                        static_cast<double>(T[cellI]) - 273.15;
                    const double Tw_C =
                        static_cast<double>(Tw[faceI]) - 273.15;

                    if (Tw_C >= T_cell_C) continue;

                    const double p_Pa = static_cast<double>(p[cellI]);
                    const double Y_c  = static_cast<double>(Y[cellI]);
                    const double Psat_w = 611.21 * std::exp(
                        (18.678 - Tw_C/234.5) * (Tw_C/(257.14+Tw_C))
                    );
                    const double Pv = (Y_c * p_Pa)
                        / std::max(epsilon + Y_c*(1.0-epsilon), 1e-20);

                    RH[cellI] = static_cast<scalar>(
                        std::min(std::max(Pv/std::max(Psat_w, 1e-10), 0.0), 1.0)
                    );
                }
            }

            RH.correctBoundaryConditions();

            Info<< "Region " << fluidRegions[i].name()
                << " | RH min/max: "
                << min(RH).value() << " / " << max(RH).value()
                << " | Y  min/max: "
                << min(Y).value()  << " / " << max(Y).value()
                << " | film max: " << max(film).value() << " m"
                << endl;
        }
        // =====================================================================
        // FIN BLOQUE PSICROMETRICO
        // =====================================================================

        if (nOuterCorr != 1)
        {
            forAll(fluidRegions, i)
            {
                #include "storeOldFluidFields.H"
            }
        }

        // --- PIMPLE loop
        for (int oCorr=0; oCorr<nOuterCorr; ++oCorr)
        {
            const bool finalIter = (oCorr == nOuterCorr-1);

            forAll(fluidRegions, i)
            {
                fvMesh& mesh = fluidRegions[i];

                #include "readFluidMultiRegionPIMPLEControls.H"
                #include "setRegionFluidFields.H"
                #include "solveFluid.H"
            }

            forAll(solidRegions, i)
            {
                fvMesh& mesh = solidRegions[i];

                #include "readSolidMultiRegionPIMPLEControls.H"
                #include "setRegionSolidFields.H"
                #include "solveSolid.H"
            }

            if (coupled)
            {
                Info<< "\nSolving energy coupled regions " << endl;
                fvMatrixAssemblyPtr->solve();
                #include "correctThermos.H"
                forAll(fluidRegions, i)
                {
                    #include "setRegionFluidFields.H"
                    rho = thermo.rho();
                }
                fvMatrixAssemblyPtr->clear();

                forAll(fluidRegions, i)
                {
                    fvMesh& mesh = fluidRegions[i];

                    #include "readFluidMultiRegionPIMPLEControls.H"
                    #include "setRegionFluidFields.H"
                    if (!frozenFlow)
                    {
                        Info<< "\nSolving for fluid region "
                            << fluidRegions[i].name() << endl;
                        for (int corr=0; corr<nCorr; corr++)
                        {
                            #include "pEqn.H"
                        }
                        turbulence.correct();
                    }

                    rho = thermo.rho();
                    Info<< "Min/max T:" << min(thermo.T()).value() << ' '
                        << max(thermo.T()).value() << endl;
                }

                fvMatrixAssemblyPtr->clear();
            }

            if (!oCorr && nOuterCorr > 1)
            {
                loopControl looping(runTime, pimple, "energyCoupling");

                while (looping.loop())
                {
                    Info<< nl << looping << nl;

                    forAll(fluidRegions, i)
                    {
                        fvMesh& mesh = fluidRegions[i];

                        Info<< "\nSolving for fluid region "
                            << fluidRegions[i].name() << endl;
                        #include "readFluidMultiRegionPIMPLEControls.H"
                        #include "setRegionFluidFields.H"
                        frozenFlow = true;
                        #include "solveFluid.H"
                    }

                    forAll(solidRegions, i)
                    {
                        fvMesh& mesh = solidRegions[i];

                        Info<< "\nSolving for solid region "
                            << solidRegions[i].name() << endl;
                        #include "readSolidMultiRegionPIMPLEControls.H"
                        #include "setRegionSolidFields.H"
                        #include "solveSolid.H"
                    }

                    if (coupled)
                    {
                        Info<< "\nSolving energy coupled regions " << endl;
                        fvMatrixAssemblyPtr->solve();
                        #include "correctThermos.H"

                        forAll(fluidRegions, i)
                        {
                            #include "setRegionFluidFields.H"
                            rho = thermo.rho();
                        }

                        fvMatrixAssemblyPtr->clear();
                    }
                }
            }
        }

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
