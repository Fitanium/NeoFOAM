// SPDX-License-Identifier: MIT
// SPDX-FileCopyrightText: 2023-2024 NeoFOAM authors
#define CATCH_CONFIG_RUNNER // Define this before including catch.hpp to create
                            // a custom main

#include "catch2_common.hpp"

#include "common.hpp"

namespace dsl = NeoFOAM::dsl;

// TEST_CASE("TemporalOperator")
TEMPLATE_TEST_CASE("TemporalOperator", "[template]", NeoFOAM::scalar, NeoFOAM::Vector)
{
    auto [execName, exec] = GENERATE(allAvailableExecutor());

    auto mesh = NeoFOAM::createSingleCellMesh(exec);
    auto sp = NeoFOAM::finiteVolume::cellCentred::SparsityPattern {mesh};

    SECTION("Operator creation on " + execName)
    {
        NeoFOAM::Field<TestType> fA(exec, 1, 2.0 * NeoFOAM::one<TestType>());
        NeoFOAM::BoundaryFields<TestType> bf(exec, mesh.boundaryMesh().offset());

        std::vector<fvcc::VolumeBoundary<TestType>> bcs {};
        auto vf = fvcc::VolumeField<TestType>(exec, "vf", mesh, fA, bf, bcs);
        dsl::TemporalOperator<TestType> b = TemporalDummy<TestType>(vf);

        REQUIRE(b.getName() == "TemporalDummy");
        REQUIRE(b.getType() == dsl::Operator::Type::Explicit);
    }

    SECTION("Supports Coefficients Explicit " + execName)
    {
        std::vector<fvcc::VolumeBoundary<TestType>> bcs {};
        NeoFOAM::scalar t = 0.0;
        NeoFOAM::scalar dt = 0.1;

        NeoFOAM::Field<TestType> fA(exec, 1, 2.0 * NeoFOAM::one<TestType>());
        NeoFOAM::Field<NeoFOAM::scalar> scaleField(exec, 1, 2.0);
        NeoFOAM::BoundaryFields<TestType> bf(exec, mesh.boundaryMesh().offset());
        auto vf = fvcc::VolumeField<TestType>(exec, "vf", mesh, fA, bf, bcs);

        dsl::TemporalOperator<TestType> c =
            2 * dsl::TemporalOperator<TestType>(TemporalDummy<TestType>(vf));
        dsl::TemporalOperator<TestType> d =
            scaleField * dsl::TemporalOperator<TestType>(TemporalDummy<TestType>(vf));
        dsl::TemporalOperator<TestType> e =
            Coeff(-3, scaleField) * dsl::TemporalOperator<TestType>(TemporalDummy<TestType>(vf));

        [[maybe_unused]] auto coeffC = c.getCoefficient();
        [[maybe_unused]] auto coeffD = d.getCoefficient();
        [[maybe_unused]] auto coeffE = e.getCoefficient();

        NeoFOAM::Field<TestType> source(exec, 1, 2.0 * NeoFOAM::one<TestType>());
        c.explicitOperation(source, t, dt);

        // 2 += 2 * 2
        auto hostSourceC = source.copyToHost();
        REQUIRE(hostSourceC.view()[0] == 6.0 * NeoFOAM::one<TestType>());

        // 6 += 2 * 2
        d.explicitOperation(source, t, dt);
        auto hostSourceD = source.copyToHost();
        REQUIRE(hostSourceD.view()[0] == 10.0 * NeoFOAM::one<TestType>());

        // 10 += - 6 * 2
        e.explicitOperation(source, t, dt);
        auto hostSourceE = source.copyToHost();
        REQUIRE(hostSourceE.view()[0] == -2.0 * NeoFOAM::one<TestType>());
    }

    SECTION("Implicit Operations " + execName)
    {
        NeoFOAM::Field<TestType> fA(exec, 1, 2.0 * NeoFOAM::one<TestType>());
        NeoFOAM::BoundaryFields<TestType> bf(exec, mesh.boundaryMesh().offset());

        std::vector<fvcc::VolumeBoundary<TestType>> bcs {};
        auto vf = fvcc::VolumeField<TestType>(exec, "vf", mesh, fA, bf, bcs);
        dsl::TemporalOperator<TestType> b = TemporalDummy<TestType>(vf, Operator::Type::Implicit);

        REQUIRE(b.getName() == "TemporalDummy");
        REQUIRE(b.getType() == Operator::Type::Implicit);
    }

    auto ls = NeoFOAM::la::createEmptyLinearSystem<
        TestType,
        NeoFOAM::localIdx,
        NeoFOAM::finiteVolume::cellCentred::SparsityPattern>(sp);

    SECTION("Supports Coefficients Implicit " + execName)
    {
        std::vector<fvcc::VolumeBoundary<TestType>> bcs {};
        NeoFOAM::scalar t = 0.0;
        NeoFOAM::scalar dt = 0.1;

        NeoFOAM::Field<TestType> fA(exec, 1, 2.0 * NeoFOAM::one<TestType>());
        NeoFOAM::Field<NeoFOAM::scalar> scaleField(exec, 1, 2.0);
        NeoFOAM::BoundaryFields<TestType> bf(exec, mesh.boundaryMesh().offset());
        auto vf = fvcc::VolumeField<TestType>(exec, "vf", mesh, fA, bf, bcs);

        auto c = 2 * dsl::TemporalOperator<TestType>(TemporalDummy(vf, Operator::Type::Implicit));
        auto d = scaleField
               * dsl::TemporalOperator<TestType>(TemporalDummy(vf, Operator::Type::Implicit));
        auto e = Coeff(-3, scaleField)
               * dsl::TemporalOperator<TestType>(TemporalDummy(vf, Operator::Type::Implicit));

        [[maybe_unused]] auto coeffC = c.getCoefficient();
        [[maybe_unused]] auto coeffD = d.getCoefficient();
        [[maybe_unused]] auto coeffE = e.getCoefficient();

        Field source(exec, 1, 2.0);
        c.implicitOperation(ls, t, dt);

        // c = 2 * 2
        auto hostRhsC = ls.rhs().copyToHost();
        REQUIRE(hostRhsC.view()[0] == 4.0 * NeoFOAM::one<TestType>());
        auto hostLsC = ls.copyToHost();
        REQUIRE(hostLsC.matrix().values().view()[0] == 4.0 * NeoFOAM::one<TestType>());


        // // d= 2 * 2
        ls.reset();
        d.implicitOperation(ls, t, dt);
        auto hostRhsD = ls.rhs().copyToHost();
        REQUIRE(hostRhsD.view()[0] == 4.0 * NeoFOAM::one<TestType>());
        auto hostLsD = ls.copyToHost();
        REQUIRE(hostLsD.matrix().values().view()[0] == 4.0 * NeoFOAM::one<TestType>());


        // e = - -3 * 2 * 2 = -12
        ls.reset();
        e.implicitOperation(ls, t, dt);
        auto hostRhsE = ls.rhs().copyToHost();
        REQUIRE(hostRhsE.view()[0] == -12.0 * NeoFOAM::one<TestType>());
        auto hostLsE = ls.copyToHost();
        REQUIRE(hostLsE.matrix().values().view()[0] == -12.0 * NeoFOAM::one<TestType>());
    }
}
