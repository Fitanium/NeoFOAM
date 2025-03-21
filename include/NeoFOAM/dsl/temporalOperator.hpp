// SPDX-License-Identifier: MIT
// SPDX-FileCopyrightText: 2023-2025 NeoFOAM authors
#pragma once

#include <memory>
#include <concepts>

#include "NeoFOAM/core/primitives/scalar.hpp"
#include "NeoFOAM/fields/field.hpp"
#include "NeoFOAM/linearAlgebra/linearSystem.hpp"
#include "NeoFOAM/core/input.hpp"
#include "NeoFOAM/dsl/coeff.hpp"
#include "NeoFOAM/dsl/operator.hpp"

namespace NeoFOAM::dsl
{

template<typename T>
concept HasTemporalExplicitOperator = requires(T t) {
    {
        t.explicitOperation(
            std::declval<Field<typename T::FieldValueType>&>(),
            std::declval<NeoFOAM::scalar>(),
            std::declval<NeoFOAM::scalar>()
        )
    } -> std::same_as<void>; // Adjust return type and arguments as needed
};

template<typename T>
concept HasTemporalImplicitOperator = requires(T t) {
    {
        t.implicitOperation(
            std::declval<la::LinearSystem<typename T::FieldValueType, localIdx>&>(),
            std::declval<NeoFOAM::scalar>(),
            std::declval<NeoFOAM::scalar>()
        )
    } -> std::same_as<void>; // Adjust return type and arguments as needed
};

template<typename T>
concept HasTemporalOperator = HasTemporalExplicitOperator<T> || HasTemporalImplicitOperator<T>;

/* @class TemporalOperator
 * @brief A class to represent a TemporalOperator in NeoFOAMs DSL
 *
 * The design here is based on the type erasure design pattern
 * see https://www.youtube.com/watch?v=4eeESJQk-mw
 *
 * Motivation for using type erasure is that concrete implementation
 * of TemporalOperator e.g Divergence, Laplacian, etc can be stored in a vector of
 * TemporalOperator
 *
 * @ingroup dsl
 */
template<typename ValueType>
class TemporalOperator
{
public:

    using FieldValueType = ValueType;

    template<HasTemporalOperator T>
    TemporalOperator(T cls) : model_(std::make_unique<TemporalOperatorModel<T>>(std::move(cls)))
    {}

    TemporalOperator(const TemporalOperator& eqnOperator) : model_ {eqnOperator.model_->clone()} {}

    TemporalOperator(TemporalOperator&& eqnOperator) : model_ {std::move(eqnOperator.model_)} {}

    void explicitOperation(Field<ValueType>& source, scalar t, scalar dt) const
    {
        model_->explicitOperation(source, t, dt);
    }

    void implicitOperation(la::LinearSystem<ValueType, localIdx>& ls, scalar t, scalar dt)
    {
        model_->implicitOperation(ls, t, dt);
    }

    la::LinearSystem<ValueType, localIdx> createEmptyLinearSystem() const
    {
        return model_->createEmptyLinearSystem();
    }

    /* returns the fundamental type of an operator, ie explicit, implicit */
    Operator::Type getType() const { return model_->getType(); }

    std::string getName() const { return model_->getName(); }

    Coeff& getCoefficient() { return model_->getCoefficient(); }

    Coeff getCoefficient() const { return model_->getCoefficient(); }

    /* @brief Given an input this function reads required coeffs */
    void build(const Input& input) { model_->build(input); }

    /* @brief Get the executor */
    const Executor& exec() const { return model_->exec(); }


private:

    /* @brief Base class defining the concept of a term. This effectively
     * defines what functions need to be implemented by a concrete Operator implementation
     * */
    struct TemporalOperatorConcept
    {
        virtual ~TemporalOperatorConcept() = default;

        virtual void explicitOperation(Field<ValueType>& source, scalar t, scalar dt) = 0;

        virtual void
        implicitOperation(la::LinearSystem<ValueType, localIdx>& ls, scalar t, scalar dt) = 0;

        virtual la::LinearSystem<ValueType, localIdx> createEmptyLinearSystem() const = 0;

        /* @brief Given an input this function reads required coeffs */
        virtual void build(const Input& input) = 0;

        /* returns the name of the operator */
        virtual std::string getName() const = 0;

        /* returns the fundamental type of an operator, ie explicit, implicit, temporal */
        virtual Operator::Type getType() const = 0;

        /* @brief get the associated coefficient for this term */
        virtual Coeff& getCoefficient() = 0;

        /* @brief get the associated coefficient for this term */
        virtual Coeff getCoefficient() const = 0;

        /* @brief Get the executor */
        virtual const Executor& exec() const = 0;

        // The Prototype Design Pattern
        virtual std::unique_ptr<TemporalOperatorConcept> clone() const = 0;
    };

    // Templated derived class to implement the type-specific behavior
    template<typename ConcreteTemporalOperatorType>
    struct TemporalOperatorModel : TemporalOperatorConcept
    {
        /* @brief build with concrete TemporalOperator */
        TemporalOperatorModel(ConcreteTemporalOperatorType concreteOp)
            : concreteOp_(std::move(concreteOp))
        {}

        /* returns the name of the operator */
        std::string getName() const override { return concreteOp_.getName(); }

        virtual void explicitOperation(Field<ValueType>& source, scalar t, scalar dt) override
        {
            if constexpr (HasTemporalExplicitOperator<ConcreteTemporalOperatorType>)
            {
                concreteOp_.explicitOperation(source, t, dt);
            }
        }

        virtual void
        implicitOperation(la::LinearSystem<ValueType, localIdx>& ls, scalar t, scalar dt) override
        {
            if constexpr (HasTemporalImplicitOperator<ConcreteTemporalOperatorType>)
            {
                concreteOp_.implicitOperation(ls, t, dt);
            }
        }

        virtual la::LinearSystem<ValueType, localIdx> createEmptyLinearSystem() const override
        {
            // if constexpr (HasTemporalImplicitOperator<ConcreteTemporalOperatorType>)
            // {
            return concreteOp_.createEmptyLinearSystem();
            // }
            throw std::runtime_error("Implicit operation not implemented");
            // only need to avoid compiler warning about missing return statement
            // this code path should never be reached as we call implicitOperation on an explicit
            // operator
            NeoFOAM::Field<ValueType> values(exec(), 1, zero<ValueType>());
            NeoFOAM::Field<NeoFOAM::localIdx> colIdx(exec(), 1, 0);
            NeoFOAM::Field<NeoFOAM::localIdx> rowPtrs(exec(), 2, 0);
            NeoFOAM::la::CSRMatrix<ValueType, NeoFOAM::localIdx> csrMatrix(values, colIdx, rowPtrs);

            NeoFOAM::Field<ValueType> rhs(exec(), 1, zero<ValueType>());
            return la::LinearSystem<ValueType, localIdx>(csrMatrix, rhs, "custom");
        }

        /* @brief Given an input this function reads required coeffs */
        virtual void build(const Input& input) override { concreteOp_.build(input); }

        /* returns the fundamental type of an operator, ie explicit, implicit, temporal */
        Operator::Type getType() const override { return concreteOp_.getType(); }

        /* @brief Get the executor */
        const Executor& exec() const override { return concreteOp_.exec(); }

        /* @brief get the associated coefficient for this term */
        virtual Coeff& getCoefficient() override { return concreteOp_.getCoefficient(); }

        /* @brief get the associated coefficient for this term */
        virtual Coeff getCoefficient() const override { return concreteOp_.getCoefficient(); }

        // The Prototype Design Pattern
        std::unique_ptr<TemporalOperatorConcept> clone() const override
        {
            return std::make_unique<TemporalOperatorModel>(*this);
        }

        ConcreteTemporalOperatorType concreteOp_;
    };

    std::unique_ptr<TemporalOperatorConcept> model_;
};


template<typename ValueType>
TemporalOperator<ValueType> operator*(scalar scalarCoeff, TemporalOperator<ValueType> rhs)
{
    TemporalOperator<ValueType> result = rhs;
    result.getCoefficient() *= scalarCoeff;
    return result;
}

template<typename ValueType>
TemporalOperator<ValueType>
operator*(const Field<scalar>& coeffField, TemporalOperator<ValueType> rhs)
{
    TemporalOperator<ValueType> result = rhs;
    result.getCoefficient() *= Coeff(coeffField);
    return result;
}

template<typename ValueType>
TemporalOperator<ValueType> operator*(const Coeff& coeff, TemporalOperator<ValueType> rhs)
{
    TemporalOperator<ValueType> result = rhs;
    result.getCoefficient() *= coeff;
    return result;
}


} // namespace NeoFOAM::dsl
