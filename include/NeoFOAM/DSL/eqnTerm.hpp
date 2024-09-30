// SPDX-License-Identifier: MIT
// SPDX-FileCopyrightText: 2023-2024 NeoFOAM authors
#pragma once

#include <iostream>
#include <memory>
#include <type_traits>
#include <utility>
#include <concepts>

#include "NeoFOAM/core/primitives/scalar.hpp"
#include "NeoFOAM/fields/field.hpp"
#include "NeoFOAM/core/parallelAlgorithms.hpp"
#include "NeoFOAM/finiteVolume/cellCentred/fields/volumeField.hpp"
#include "NeoFOAM/core/inputs.hpp"


namespace fvcc = NeoFOAM::finiteVolume::cellCentred;

namespace NeoFOAM::DSL
{

template<typename T>
concept HasTemporalTerm = requires(T t) {
    {
        t.temporalOperation(std::declval<Field<scalar>&>(), std::declval<scalar>())
    } -> std::same_as<void>; // Adjust return type and arguments as needed
};

template<typename T>
concept HasExplicitTerm = requires(T t) {
    {
        t.explicitOperation(std::declval<Field<scalar>&>())
    } -> std::same_as<void>; // Adjust return type and arguments as needed
};


/* @class EqnTermMixin
 * @brief A mixin class to represent a term in NeoFOAMs DSL
 *
 * @ingroup DSL
 */
template<typename ValueType>
class EqnTermMixin
{
public:

    using EqnTermValueType = ValueType;

    EqnTermMixin(bool isEvaluated) : field_(), scaleField_(1.0), termEvaluated(isEvaluated) {};

    virtual ~EqnTermMixin() = default;

    void setField(std::shared_ptr<Field<scalar>> field)
    {
        field_ = field;
        scaleField() = ScalingField<ValueType>(field_->span(), scaleField().value);
    }

    ScalingField<ValueType>& scaleField() { return scaleField_; }

    ScalingField<ValueType> scaleField() const { return scaleField_; }

    bool evaluated() { return termEvaluated; }

    std::shared_ptr<Field<scalar>>& field() { return field_; }

protected:

    std::shared_ptr<Field<scalar>> field_;

    ScalingField<ValueType> scaleField_;

    bool termEvaluated;
};


/* @class EqnTerm
 * @brief A class to represent a term in NeoFOAMs DSL
 *
 * @ingroup DSL
 */
template<typename ValueType>
class EqnTerm
{
public:

    enum class Type
    {
        Temporal,
        Implicit,
        Explicit
    };

    // expise valuetype to other templates
    using EqnTermValueType = ValueType;

    template<typename T>
    EqnTerm(T cls) : model_(std::make_unique<Model<T>>(std::move(cls)))
    {}

    EqnTerm(const EqnTerm& eqnTerm) : model_ {eqnTerm.model_->clone()} {}

    EqnTerm(EqnTerm&& eqnTerm) : model_ {std::move(eqnTerm.model_)} {}

    EqnTerm& operator=(const EqnTerm& eqnTerm)
    {
        model_ = eqnTerm.model_->clone();
        return *this;
    }

    EqnTerm& operator=(EqnTerm&& eqnTerm)
    {
        model_ = std::move(eqnTerm.model_);
        return *this;
    }

    std::string display() const { return model_->display(); }

    void explicitOperation(Field<ValueType>& source) { model_->explicitOperation(source); }

    void temporalOperation(Field<ValueType>& field) { model_->temporalOperation(field); }

    EqnTerm::Type getType() const { return model_->getType(); }

    void setField(std::shared_ptr<Field<scalar>> field) { model_->setField(field); }

    ScalingField<ValueType>& scaleField() { return model_->scaleField(); }

    ScalingField<ValueType> scaleField() const { return model_->scaleField(); }

    bool evaluated() { return model_->evaluated(); }

    void build(const Input& input) { model_->build(input); }

    const Executor& exec() const { return model_->exec(); }

    std::size_t nCells() const { return model_->nCells(); }

    fvcc::VolumeField<scalar>* volumeField() { return model_->volumeField(); }


private:

    // Base class to hold the type-erased value and the display function
    struct Concept
    {
        virtual ~Concept() = default;
        virtual std::string display() const = 0;
        virtual void explicitOperation(Field<scalar>& source) = 0;
        virtual void temporalOperation(Field<scalar>& field) = 0;
        virtual void build(const Input& input) = 0;

        virtual bool evaluated() = 0;

        virtual EqnTerm::Type getType() const = 0;

        virtual void setField(std::shared_ptr<Field<scalar>> field) = 0;
        virtual ScalingField<scalar>& scaleField() = 0;
        virtual ScalingField<scalar> scaleField() const = 0;

        virtual const Executor& exec() const = 0;
        virtual std::size_t nCells() const = 0;
        virtual fvcc::VolumeField<scalar>* volumeField() = 0;

        // The Prototype Design Pattern
        virtual std::unique_ptr<Concept> clone() const = 0;
    };

    // Templated derived class to implement the type-specific behavior
    template<typename T>
    struct Model : Concept
    {
        Model(T cls) : cls_(std::move(cls)) {}

        std::string display() const override { return cls_.display(); }

        virtual void explicitOperation(Field<scalar>& source) override
        {
            if constexpr (HasExplicitTerm<T>)
            {
                cls_.explicitOperation(source);
            }
        }

        virtual void temporalOperation(Field<scalar>& field) override
        {
            if constexpr (HasTemporalTerm<T>)
            {
                cls_.temporalOperation(field);
            }
        }

        virtual fvcc::VolumeField<scalar>* volumeField() override { return cls_.volumeField(); }

        virtual void build(const Input& input) override { cls_.build(input); }

        virtual bool evaluated() override { return cls_.evaluated(); }

        EqnTerm::Type getType() const override { return cls_.getType(); }

        const Executor& exec() const override { return cls_.exec(); }

        std::size_t nCells() const override { return cls_.nCells(); }

        void setField(std::shared_ptr<Field<scalar>> field) override { cls_.setField(field); }

        virtual ScalingField<scalar>& scaleField() override { return cls_.scaleField(); }

        virtual ScalingField<scalar> scaleField() const override { return cls_.scaleField(); }

        // The Prototype Design Pattern
        std::unique_ptr<Concept> clone() const override { return std::make_unique<Model>(*this); }

        T cls_;
    };

    std::unique_ptr<Concept> model_;
};


// add multiply operator to EqnTerm
template<typename EqnTermType>
auto operator*(scalar scale, const EqnTermType& lhs)
{
    using ValueType = typename EqnTermType::EqnTermValueType;
    EqnTerm<ValueType> result = lhs;
    result.scaleField() *= scale;
    return result;
}

// add multiply operator to EqnTerm
template<typename EqnTermType>
auto operator*(Field<scalar> scale, const EqnTermType& lhs)
{
    using ValueType = typename EqnTermType::EqnTermValueType;
    EqnTerm<ValueType> result = lhs;
    if (!result.scaleField().useSpan)
    {
        // allocate the scaling field to avoid deallocation
        result.setField(std::make_shared<Field<scalar>>(scale));
    }
    else
    {
        result.scaleField() *= scale;
    }
    return result;
}

template<typename ValueType, parallelForFieldKernel<ValueType> ScaleFunction>
EqnTerm<ValueType> operator*(ScaleFunction scaleFunc, const EqnTerm<ValueType>& lhs)
{
    EqnTerm<ValueType> result = lhs;
    if (!result.scaleField().useSpan)
    {
        result.setField(std::make_shared<Field<scalar>>(result.exec(), result.nCells(), 1.0));
    }
    map(result.exec(), result.scaleField().values, scaleFunc);
    return result;
}

} // namespace NeoFOAM::DSL
