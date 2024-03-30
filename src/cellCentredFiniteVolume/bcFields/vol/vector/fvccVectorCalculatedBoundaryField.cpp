#include "NeoFOAM/cellCentredFiniteVolume/bcFields/vol/vector/fvccVectorCalculatedBoundaryField.hpp"

namespace NeoFOAM
{
fvccVectorCalculatedBoundaryField::fvccVectorCalculatedBoundaryField(const unstructuredMesh& mesh, int patchID)
    : fvccBoundaryField<Vector>(mesh, patchID)
{
}

void fvccVectorCalculatedBoundaryField::correctBoundaryConditions(boundaryFields<Vector>& bfield, const Field<Vector>& internalField)
{

}


}