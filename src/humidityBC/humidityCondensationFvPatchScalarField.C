/*---------------------------------------------------------------------------*\
    PASO 2: implementacion con logica termodinamica.

    En cada cara del patch:
      - Lee T_wall y p en la celda
      - Calcula Y_eq = e*Psat(T_wall) / (p - (1-e)*Psat)   (Buck + Eq.2.4 TFM)
      - Si Y_celda > Y_eq -> condensa  -> Y_wall = Y_eq (Dirichlet)
      - Si Y_celda <= Y_eq -> no condensa -> zeroGradient
\*---------------------------------------------------------------------------*/

#include "humidityCondensationFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

namespace Foam
{

// * * * * * * * * * * * * * * * Constructores * * * * * * * * * * * * * * * //

humidityCondensationFvPatchScalarField::humidityCondensationFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    TName_("T"),
    pName_("p"),
    epsilon_(0.622)
{
    this->refValue()      = scalarField(p.size(), 0.0);
    this->refGrad()       = scalarField(p.size(), 0.0);
    this->valueFraction() = scalarField(p.size(), 0.0);
}


humidityCondensationFvPatchScalarField::humidityCondensationFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    TName_(dict.getOrDefault<word>("T", "T")),
    pName_(dict.getOrDefault<word>("p", "p")),
    epsilon_(dict.getOrDefault<scalar>("epsilon", 0.622))
{
    if (dict.found("refValue"))
    {
        this->refValue() = scalarField("refValue", dict, p.size());
    }
    else
    {
        this->refValue() = scalarField(p.size(), 0.0);
    }

    if (dict.found("refGradient"))
    {
        this->refGrad() = scalarField("refGradient", dict, p.size());
    }
    else
    {
        this->refGrad() = scalarField(p.size(), 0.0);
    }

    if (dict.found("valueFraction"))
    {
        this->valueFraction() = scalarField("valueFraction", dict, p.size());
    }
    else
    {
        this->valueFraction() = scalarField(p.size(), 0.0);
    }

    if (dict.found("value"))
    {
        fvPatchScalarField::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        this->evaluate();
    }
}


humidityCondensationFvPatchScalarField::humidityCondensationFvPatchScalarField
(
    const humidityCondensationFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    TName_(ptf.TName_),
    pName_(ptf.pName_),
    epsilon_(ptf.epsilon_)
{}


humidityCondensationFvPatchScalarField::humidityCondensationFvPatchScalarField
(
    const humidityCondensationFvPatchScalarField& ptf
)
:
    mixedFvPatchScalarField(ptf),
    TName_(ptf.TName_),
    pName_(ptf.pName_),
    epsilon_(ptf.epsilon_)
{}


humidityCondensationFvPatchScalarField::humidityCondensationFvPatchScalarField
(
    const humidityCondensationFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    TName_(ptf.TName_),
    pName_(ptf.pName_),
    epsilon_(ptf.epsilon_)
{}


// * * * * * * * * * * * * * * * Metodos miembro * * * * * * * * * * * * * * //

void humidityCondensationFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // -----------------------------------------------------------------------
    // Acceso a campos del registry
    // -----------------------------------------------------------------------
    // T en este patch (valor en cara = T_pared, gracias al CHT acoplado)
    const fvPatchField<scalar>& Tp =
        patch().lookupPatchField<volScalarField, scalar>(TName_);
    // p en este patch (gradiente bajo cerca de pared, p_cara ~ p_celda)
    const fvPatchField<scalar>& pp =
        patch().lookupPatchField<volScalarField, scalar>(pName_);

    // Y en la celda interior adyacente al patch (valor de la celda, no de la cara)
    const scalarField Yc(this->patchInternalField());

    // Referencias modificables a los campos de la BC mixta
    scalarField& vRef  = this->refValue();
    scalarField& vGrad = this->refGrad();
    scalarField& vFrac = this->valueFraction();

    // -----------------------------------------------------------------------
    // Bucle cara a cara: termodinamica psicrometrica
    // -----------------------------------------------------------------------
    forAll(Tp, faceI)
    {
        const scalar Tw_C = Tp[faceI] - 273.15;   // T pared en Celsius
        const scalar p_Pa = pp[faceI];            // p local
        const scalar Y_c  = Yc[faceI];            // Y en celda adyacente

        // Buck (1981): presion de saturacion a T_pared [Pa]
        const scalar Psat_w = 611.21 * Foam::exp(
            (18.678 - Tw_C / 234.5) * (Tw_C / (257.14 + Tw_C))
        );

        // Eq. 2.4 TFM: Y_eq = e*Psat / (p - (1-e)*Psat)
        // Esto es la fraccion masica del vapor en saturacion (RH=100%) a T_pared
        const scalar Y_eq = epsilon_ * Psat_w
            / max(p_Pa - (1.0 - epsilon_) * Psat_w, scalar(1.0));

        if (Y_c > Y_eq)
        {
            // ------------------------------------------------------------
            // Sobresaturacion -> CONDENSA
            // Imponemos Dirichlet: Y_pared = Y_eq(T_pared)
            // El gradiente que aparece chupa vapor del aire => mDot natural
            // ------------------------------------------------------------
            vRef[faceI]  = Y_eq;
            vGrad[faceI] = 0.0;
            vFrac[faceI] = 1.0;
        }
        else
        {
            // ------------------------------------------------------------
            // No hay sobresaturacion -> NO CONDENSA
            // zeroGradient: la pared no aporta ni quita vapor
            // ------------------------------------------------------------
            vRef[faceI]  = 0.0;
            vGrad[faceI] = 0.0;
            vFrac[faceI] = 0.0;
        }
    }

    mixedFvPatchScalarField::updateCoeffs();
}


void humidityCondensationFvPatchScalarField::write(Ostream& os) const
{
    mixedFvPatchScalarField::write(os);
    os.writeEntry("T", TName_);
    os.writeEntry("p", pName_);
    os.writeEntry("epsilon", epsilon_);
}


// * * * * * * * * * * * * * Registro en runtime * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    humidityCondensationFvPatchScalarField
);

} // End namespace Foam
