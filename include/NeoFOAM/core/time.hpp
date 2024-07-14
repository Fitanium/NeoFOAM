// SPDX-License-Identifier: MIT
// SPDX-FileCopyrightText: 2023 NeoFOAM authors
#pragma once

#include <iostream>
#include <string>

#include "primitives/scalar.hpp"

const std::string NL = "\n";

namespace NeoFOAM
{

using word = std::string;

class ArgList
{

public:

    ArgList(int argc, char* argv[]) {};

    [[nodiscard]] bool checkRootCase() const { return true; };
};

class Time
{
public:

    const static word CONTROL_DICT_NAME;

    Time() : time_(0.) {};

    Time(const word, const ArgList) : time_(0.) {};

    [[nodiscard]] word timeName() { return std::to_string(time_); }

    [[nodiscard]] bool loop()
    {
        time_ += 1.;
        return 10.0 >= time_;
    };

    std::ostream& printExecutionTime(std::ostream& os) const { return os; };

private:

    scalar time_;
};

const word Time::CONTROL_DICT_NAME = "system/controlDict";

} // namespace NeoFOAM
