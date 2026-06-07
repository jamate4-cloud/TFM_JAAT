// Initialise fluid field pointer lists
PtrList<rhoReactionThermo> thermoFluid(fluidRegions.size());
PtrList<volScalarField> rhoFluid(fluidRegions.size());
PtrList<volVectorField> UFluid(fluidRegions.size());
PtrList<volScalarField> RHFluid(fluidRegions.size());
PtrList<volScalarField> YFluid(fluidRegions.size());          // Fracción másica vapor [-], calculada internamente
PtrList<volScalarField> filmThicknessFluid(fluidRegions.size());
PtrList<surfaceScalarField> phiFluid(fluidRegions.size());
PtrList<uniformDimensionedScalarField> hRefFluid(fluidRegions.size());
PtrList<volScalarField> ghFluid(fluidRegions.size());
PtrList<surfaceScalarField> ghfFluid(fluidRegions.size());
PtrList<compressible::turbulenceModel> turbulenceFluid(fluidRegions.size());
PtrList<CombustionModel<rhoReactionThermo>> reactionFluid(fluidRegions.size());
PtrList<volScalarField> p_rghFluid(fluidRegions.size());
PtrList<radiation::radiationModel> radiation(fluidRegions.size());
PtrList<volScalarField> KFluid(fluidRegions.size());
PtrList<volScalarField> dpdtFluid(fluidRegions.size());
PtrList<multivariateSurfaceInterpolationScheme<scalar>::fieldTable>
    fieldsFluid(fluidRegions.size());
PtrList<volScalarField> QdotFluid(fluidRegions.size());

PtrList<fvVectorMatrix> UEqFluid(fluidRegions.size());

List<scalar> initialMassFluid(fluidRegions.size());
List<bool> frozenFlowFluid(fluidRegions.size(), false);

PtrList<IOMRFZoneList> MRFfluid(fluidRegions.size());
PtrList<fv::options> fluidFvOptions(fluidRegions.size());

List<label> pRefCellFluid(fluidRegions.size());
List<scalar> pRefValueFluid(fluidRegions.size());

PtrList<dimensionedScalar> rhoMinFluid(fluidRegions.size());
PtrList<dimensionedScalar> rhoMaxFluid(fluidRegions.size());

PtrList<pressureControl> pressureControls(fluidRegions.size());

const uniformDimensionedVectorField& g = meshObjects::gravity::New(runTime);

// Populate fluid field pointer lists
forAll(fluidRegions, i)
{
    Info<< "*** Reading fluid mesh thermophysical properties for region "
        << fluidRegions[i].name() << nl << endl;

    Info<< "    Adding to thermoFluid\n" << endl;
    thermoFluid.set(i, rhoReactionThermo::New(fluidRegions[i]).ptr());

    Info<< "    Adding to rhoFluid\n" << endl;
    rhoFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "rho",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            thermoFluid[i].rho()
        )
    );

    Info<< "    Adding to UFluid\n" << endl;
    UFluid.set
    (
        i,
        new volVectorField
        (
            IOobject
            (
                "U",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i]
        )
    );

    Info<< "    Adding to phiFluid\n" << endl;
    phiFluid.set
    (
        i,
        new surfaceScalarField
        (
            IOobject
            (
                "phi",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            linearInterpolate(rhoFluid[i]*UFluid[i])
                & fluidRegions[i].Sf()
        )
    );

    Info<< "    Adding to hRefFluid\n" << endl;
    hRefFluid.set
    (
        i,
        new uniformDimensionedScalarField
        (
            IOobject
            (
                "hRef",
                runTime.constant(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            dimensionedScalar(word::null, dimLength, Zero)
        )
    );

    dimensionedScalar ghRef
    (
        mag(g.value()) > SMALL
      ? g & (cmptMag(g.value())/mag(g.value()))*hRefFluid[i]
      : dimensionedScalar("ghRef", g.dimensions()*dimLength, 0)
    );

    Info<< "    Adding to ghFluid\n" << endl;
    ghFluid.set
    (
        i,
        new volScalarField
        (
            "gh",
            (g & fluidRegions[i].C()) - ghRef
        )
    );

    Info<< "    Adding to ghfFluid\n" << endl;
    ghfFluid.set
    (
        i,
        new surfaceScalarField
        (
            "ghf",
            (g & fluidRegions[i].Cf()) - ghRef
        )
    );

    Info<< "    Adding to turbulenceFluid\n" << endl;
    turbulenceFluid.set
    (
        i,
        compressible::turbulenceModel::New
        (
            rhoFluid[i],
            UFluid[i],
            phiFluid[i],
            thermoFluid[i]
        ).ptr()
    );

    Info<< "    Adding to reactionFluid\n" << endl;
    reactionFluid.set
    (
        i,
        CombustionModel<rhoReactionThermo>::New
        (
            thermoFluid[i],
            turbulenceFluid[i]
        )
    );

    p_rghFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "p_rgh",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i]
        )
    );

    // Force p_rgh to be consistent with p
    p_rghFluid[i] = thermoFluid[i].p() - rhoFluid[i]*ghFluid[i];

    fluidRegions[i].setFluxRequired(p_rghFluid[i].name());

    Info<< "    Adding to radiationFluid\n" << endl;
    radiation.set
    (
        i,
        radiation::radiationModel::New(thermoFluid[i].T())
    );

    initialMassFluid[i] = fvc::domainIntegrate(rhoFluid[i]).value();

    Info<< "    Adding to KFluid\n" << endl;
    KFluid.set
    (
        i,
        new volScalarField
        (
            "K",
            0.5*magSqr(UFluid[i])
        )
    );

    Info<< "    Adding to dpdtFluid\n" << endl;
    dpdtFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "dpdt",
                runTime.timeName(),
                fluidRegions[i]
            ),
            fluidRegions[i],
            dimensionedScalar(thermoFluid[i].p().dimensions()/dimTime, Zero)
        )
    );

    Info<< "    Adding to fieldsFluid\n" << endl;
    fieldsFluid.set
    (
        i,
        new multivariateSurfaceInterpolationScheme<scalar>::fieldTable
    );
    forAll(thermoFluid[i].composition().Y(), j)
    {
        fieldsFluid[i].add(thermoFluid[i].composition().Y()[j]);
    }
    fieldsFluid[i].add(thermoFluid[i].he());

    Info<< "    Adding to QdotFluid\n" << endl;
    QdotFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "Qdot",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i],
            dimensionedScalar(dimEnergy/dimVolume/dimTime, Zero)
        )
    );

    const dictionary& pimpleDict =
        fluidRegions[i].solutionDict().subDict("PIMPLE");
    pimpleDict.readIfPresent("frozenFlow", frozenFlowFluid[i]);

    rhoMaxFluid.set
    (
        i,
        new dimensionedScalar("rhoMax", dimDensity, GREAT, pimpleDict)
    );

    rhoMinFluid.set
    (
        i,
        new dimensionedScalar("rhoMin", dimDensity, Zero, pimpleDict)
    );

    pressureControls.set
    (
        i,
        new pressureControl(thermoFluid[i].p(), rhoFluid[i], pimpleDict, false)
    );

    Info<< "    Adding MRF\n" << endl;
    MRFfluid.set
    (
        i,
        new IOMRFZoneList(fluidRegions[i])
    );

    Info<< "    Adding fvOptions\n" << endl;
    fluidFvOptions.set
    (
        i,
        new fv::options(fluidRegions[i])
    );

    UEqFluid.set
    (
        i,
        new fvVectorMatrix(UFluid[i], dimForce)
    );

    turbulenceFluid[i].validate();

    pRefCellFluid[i] = -1;
    pRefValueFluid[i] = 0.0;

    if (p_rghFluid[i].needReference())
    {
        setRefCell
        (
            thermoFluid[i].p(),
            p_rghFluid[i],
            pimpleDict,
            pRefCellFluid[i],
            pRefValueFluid[i]
        );
    }

    // --- filmThickness: leído de carpeta 0 (valor inicial = 0 en todo el dominio)
    Info<< "    Adding to filmThicknessFluid\n" << endl;
    filmThicknessFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "filmThickness",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i]
        )
    );

    // --- RH: leído de carpeta 0 como condición inicial (e.g. 0.8 = 80%)
    Info<< "    Adding to RHFluid\n" << endl;
    RHFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "RH",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i]
        )
    );

    // --- Y_vapor: fracción másica de vapor [kg_vapor/kg_mezcla], dimensionless [-]
    //     NO se lee de carpeta 0. Se calcula aquí desde T, p y RH usando:
    //       Psat = Buck equation (Eq. 2.2 TFM)
    //       Pv   = RH * Psat
    //       Y    = 0.622 * Pv / (p - 0.378 * Pv)   (Eq. 2.4 TFM, epsilon=0.622)
    //
    //     BCs:
    //       - AIR_INT_to_LENS -> humidityCondensation (BC custom):
    //           impone Y = Y_eq(T_pared) cuando hay sobresaturacion,
    //           permitiendo deplecion termodinamicamente correcta.
    //       - resto           -> zeroGradient
    Info<< "    Adding to YFluid (derived from T, p, RH via psychrometric equations)\n" << endl;
    {
        const volScalarField& T_init  = thermoFluid[i].T();
        const volScalarField& p_init  = thermoFluid[i].p();
        const volScalarField& RH_init = RHFluid[i];

        // Construir lista de tipos de BC para cada patch
        wordList YPatchTypes
        (
            fluidRegions[i].boundary().size(),
            zeroGradientFvPatchScalarField::typeName
        );

        forAll(fluidRegions[i].boundary(), patchI)
        {
            const word& patchName =
                fluidRegions[i].boundary()[patchI].name();
            if (patchName == "AIR_INT_to_LENS")
            {
                YPatchTypes[patchI] = "humidityCondensation";
                Info<< "    -> Patch '" << patchName
                    << "' set to humidityCondensation BC for Y_vapor"
                    << endl;
            }
        }

        YFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "Y_vapor",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                fluidRegions[i],
                dimensionedScalar("Y_vapor", dimless, 0.0),
                YPatchTypes
            )
        );

        volScalarField& Y = YFluid[i];

        // Inicialización celda a celda — carta psicrométrica
        forAll(Y, cellI)
        {
            double T_C  = static_cast<double>(T_init[cellI]) - 273.15;
            double p_Pa = static_cast<double>(p_init[cellI]);
            double rh   = static_cast<double>(RH_init[cellI]);

            // Buck (1981): Psat [Pa]
            double Psat = 611.21 * std::exp(
                (18.678 - T_C / 234.5) * (T_C / (257.14 + T_C))
            );

            double Pv = rh * Psat;

            // Eq. 2.4 TFM: Y = 0.622*Pv / (p - 0.378*Pv)
            Y[cellI] = static_cast<scalar>(
                0.622 * Pv / std::max(p_Pa - 0.378 * Pv, 1.0)
            );
        }
        Y.correctBoundaryConditions();

        Info<< "    Y_vapor initialised. Min/Max: "
            << min(Y).value() << " / " << max(Y).value() << endl;
    }

    Info<< "    All fields initialised for region "
        << fluidRegions[i].name() << nl << endl;
}
