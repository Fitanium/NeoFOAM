// SPDX-License-Identifier: MIT
// SPDX-FileCopyrightText: 2023 NeoFOAM authors

#include "NeoFOAM/dsl/coeff.hpp"

namespace NeoFOAM::dsl
{

Coeff::Coeff() : coeff_(1.0), view_(), hasView_(false) {}

Coeff::Coeff(scalar value) : coeff_(value), view_(), hasView_(false) {}

Coeff::Coeff(scalar coeff, const Field<scalar>& field)
    : coeff_(coeff), view_(field.view()), hasView_(true)
{}

Coeff::Coeff(const Field<scalar>& field) : coeff_(1.0), view_(field.view()), hasView_(true) {}

bool Coeff::hasSpan() { return hasView_; }

View<const scalar> Coeff::span() { return view_; }

Coeff& Coeff::operator*=(scalar rhs)
{
    coeff_ *= rhs;
    return *this;
}

Coeff& Coeff::operator*=(const Coeff& rhs)
{
    if (hasView_ && rhs.hasView_)
    {
        NF_ERROR_EXIT("Not implemented");
    }

    if (!hasView_ && rhs.hasView_)
    {
        // Take over the span
        view_ = rhs.view_;
        hasView_ = true;
    }

    return this->operator*=(rhs.coeff_);
}

namespace detail
{
void toField(Coeff& coeff, Field<scalar>& rhs)
{
    if (coeff.hasSpan())
    {
        rhs.resize(coeff.span().size());
        fill(rhs, 1.0);
        auto rhsView = rhs.view();
        // otherwise we are unable to capture values in the lambda
        parallelFor(
            rhs.exec(), rhs.range(), KOKKOS_LAMBDA(const size_t i) { rhsView[i] *= coeff[i]; }
        );
    }
    else
    {
        rhs.resize(1);
        fill(rhs, coeff[0]);
    }
}
} // namespace detail


} // namespace NeoFOAM::dsl
