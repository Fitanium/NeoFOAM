// SPDX-License-Identifier: MIT
// SPDX-FileCopyrightText: 2023 NeoFOAM authors

#pragma once

#include <functional>

#include <Kokkos_Core.hpp>

#include "NeoFOAM/core/executor/executor.hpp"
#include "NeoFOAM/core/input.hpp"
#include "NeoFOAM/core/runtimeSelectionFactory.hpp"
#include "NeoFOAM/mesh/unstructured.hpp"
#include "NeoFOAM/finiteVolume/cellCentred/fields/surfaceField.hpp"
#include "NeoFOAM/finiteVolume/cellCentred/fields/volumeField.hpp"
#include "NeoFOAM/finiteVolume/cellCentred/boundary.hpp"

namespace NeoFOAM::finiteVolume::cellCentred
{

class FaceNormalGradientFactory :
    public NeoFOAM::RuntimeSelectionFactory<
        FaceNormalGradientFactory,
        Parameters<const Executor&, const UnstructuredMesh&, const Input&>>
{

public:

    static std::unique_ptr<FaceNormalGradientFactory>
    create(const Executor& exec, const UnstructuredMesh& uMesh, const Input& inputs)
    {
        // input is dictionary the key is "interpolation"
        std::string key =
            (std::holds_alternative<NeoFOAM::Dictionary>(inputs))
                ? std::get<NeoFOAM::Dictionary>(inputs).get<std::string>("faceNormalGradient")
                : std::get<NeoFOAM::TokenList>(inputs).next<std::string>();

        keyExistsOrError(key);
        return table().at(key)(exec, uMesh, inputs);
    }

    static std::string name() { return "FaceNormalGradientFactory"; }

    FaceNormalGradientFactory(const Executor& exec, const UnstructuredMesh& mesh)
        : exec_(exec), mesh_(mesh) {};

    virtual ~FaceNormalGradientFactory() {} // Virtual destructor

    virtual void faceNormalGrad(
        const VolumeField<scalar>& volField, SurfaceField<scalar>& surfaceField
    ) const = 0;

    virtual const SurfaceField<scalar>& deltaCoeffs() const = 0;

    // Pure virtual function for cloning
    virtual std::unique_ptr<FaceNormalGradientFactory> clone() const = 0;

protected:

    const Executor exec_;
    const UnstructuredMesh& mesh_;
};

class FaceNormalGradient
{

public:

    FaceNormalGradient(const FaceNormalGradient& faceNGrad)
        : exec_(faceNGrad.exec_), mesh_(faceNGrad.mesh_),
          faceNormalGradKernel_(faceNGrad.faceNormalGradKernel_->clone()) {};

    FaceNormalGradient(FaceNormalGradient&& faceNGrad)
        : exec_(faceNGrad.exec_), mesh_(faceNGrad.mesh_),
          faceNormalGradKernel_(std::move(faceNGrad.faceNormalGradKernel_)) {};

    FaceNormalGradient(
        const Executor& exec,
        const UnstructuredMesh& mesh,
        std::unique_ptr<FaceNormalGradientFactory> interpolationKernel
    )
        : exec_(exec), mesh_(mesh), faceNormalGradKernel_(std::move(interpolationKernel)) {};

    FaceNormalGradient(const Executor& exec, const UnstructuredMesh& mesh, const Input& input)
        : exec_(exec), mesh_(mesh),
          faceNormalGradKernel_(FaceNormalGradientFactory::create(exec, mesh, input)) {};


    void
    faceNormalGrad(const VolumeField<scalar>& volField, SurfaceField<scalar>& surfaceField) const
    {
        faceNormalGradKernel_->faceNormalGrad(volField, surfaceField);
    }

    const SurfaceField<scalar>& deltaCoeffs() const { return faceNormalGradKernel_->deltaCoeffs(); }


    SurfaceField<scalar> faceNormalGrad(const VolumeField<scalar>& volField) const
    {
        std::string nameInterpolated = "interpolated_" + volField.name;
        SurfaceField<scalar> surfaceField(
            exec_, nameInterpolated, mesh_, createCalculatedBCs<SurfaceBoundary<scalar>>(mesh_)
        );
        faceNormalGrad(volField, surfaceField);
        return surfaceField;
    }

private:

    const Executor exec_;
    const UnstructuredMesh& mesh_;
    std::unique_ptr<FaceNormalGradientFactory> faceNormalGradKernel_;
};


} // namespace NeoFOAM
