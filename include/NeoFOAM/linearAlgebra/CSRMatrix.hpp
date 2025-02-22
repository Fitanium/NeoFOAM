// SPDX-License-Identifier: MIT
// SPDX-FileCopyrightText: 2023 NeoFOAM authors
#pragma once

#include "NeoFOAM/fields/field.hpp"


namespace NeoFOAM::la
{

enum BlockStructType
{
    cell,
    component
};


template<typename ValueType, typename IndexType>
class CSRMatrix
{

public:

    CSRMatrix(
        const Field<ValueType>& values,
        const Field<IndexType>& colIdxs,
        const Field<IndexType>& rowPtrs
    )
        : values_(values), colIdxs_(colIdxs), rowPtrs_(rowPtrs)
    {
        NF_ASSERT(values.exec() == colIdxs.exec(), "Executors are not the same");
        NF_ASSERT(values.exec() == rowPtrs.exec(), "Executors are not the same");
    };

    ~CSRMatrix() = default;

    [[nodiscard]] const Executor& exec() const { return values_.exec(); }

    [[nodiscard]] IndexType nRows() const { return rowPtrs_.size() - 1; }

    [[nodiscard]] IndexType nValues() const { return values_.size(); }

    [[nodiscard]] IndexType nColIdxs() const { return colIdxs_.size(); }

    [[nodiscard]] std::span<ValueType> values() { return values_.span(); }
    [[nodiscard]] std::span<IndexType> colIdxs() { return colIdxs_.span(); }
    [[nodiscard]] std::span<IndexType> rowPtrs() { return rowPtrs_.span(); }

    [[nodiscard]] const std::span<const ValueType> values() const { return values_.span(); }
    [[nodiscard]] const std::span<const IndexType> colIdxs() const { return colIdxs_.span(); }
    [[nodiscard]] const std::span<const IndexType> rowPtrs() const { return rowPtrs_.span(); }


    [[nodiscard]] CSRMatrix<ValueType, IndexType> copyToExecutor(Executor dstExec) const
    {
        if (dstExec == values_.exec()) return *this;
        CSRMatrix<ValueType, IndexType> result(
            values_.copyToHost(), colIdxs_.copyToHost(), rowPtrs_.copyToHost()
        );
        result.type_ = type_;
        return result;
    }

    [[nodiscard]] CSRMatrix<ValueType, IndexType> copyToHost() const
    {
        return copyToExecutor(SerialExecutor());
    }

    const ValueType& entry(const IndexType i, const IndexType j) const
    {
        const auto& vals = values();
        const auto& cols = colIdxs();
        const auto& rows = rowPtrs();
        for (IndexType ic = 0; ic < (rows[i + 1] - rows[i]); ++ic)
        {
            if (cols[rows[i] + ic] == j)
            {
                return vals[rows[i] + ic];
            }
            if (cols[rows[i] + ic] > j) break;
        }
        NF_ERROR_EXIT("Matrix entry " << i << ", " << j << " has not been allocated, cannot get.");
    }

    ValueType& entry(const IndexType i, const IndexType j)
    {
        for (IndexType iColumn = 0; iColumn < (colIdxs()[i + 1] - colIdxs()[i]); ++iColumn)
        {
            if (colIdxs()[rowPtrs()[i] + iColumn] == j)
            {
                return values()[rowPtrs()[i] + iColumn];
            }

            if (colIdxs()[rowPtrs()[i] + iColumn] > j)
            {
                // Insert
                NF_ERROR_EXIT("Not implemented must allocate " << i << ", " << j << ".");
            }
        }
    }

private:

    BlockStructType type_ {cell};
    Field<ValueType> values_;
    Field<IndexType> colIdxs_;
    Field<IndexType> rowPtrs_;
};

} // namespace NeoFOAM
