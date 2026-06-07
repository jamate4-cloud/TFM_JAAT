/*---------------------------------------------------------------------------*\
    Boundary condition for water vapour mass fraction Y on cold walls.
    
    PASO 1: skeleton only -- behaves as zeroGradient.
            Logic will be added in PASO 2.
    
    Inherits from mixedFvPatchScalarField:
      - refValue       : Y_eq(T_wall)        when condensing
      - refGrad        : 0                   when not condensing
      - valueFraction  : 1 condensing, 0 otherwise
\*---------------------------------------------------------------------------*/

#ifndef humidityCondensationFvPatchScalarField_H
#define humidityCondensationFvPatchScalarField_H

#include "mixedFvPatchFields.H"

namespace Foam
{

class humidityCondensationFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    // Private data
    
        //- Name of temperature field (default "T")
        word TName_;
        
        //- Name of pressure field (default "p")
        word pName_;
        
        //- Mw_water / Mw_air ratio (default 0.622)
        scalar epsilon_;


public:

    //- Runtime type information
    TypeName("humidityCondensation");


    // Constructors

        //- Construct from patch and internal field
        humidityCondensationFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        humidityCondensationFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping from another instance onto a new patch
        humidityCondensationFvPatchScalarField
        (
            const humidityCondensationFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy construct
        humidityCondensationFvPatchScalarField
        (
            const humidityCondensationFvPatchScalarField&
        );

        //- Construct as copy with new internal field
        humidityCondensationFvPatchScalarField
        (
            const humidityCondensationFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new humidityCondensationFvPatchScalarField(*this)
            );
        }

        //- Clone with new internal field
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new humidityCondensationFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        //- Update the coefficients (computes Y_eq, sets BC type)
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};

} // End namespace Foam

#endif
