#ifndef EDIT_DISTANCE_REF_HPP
#define EDIT_DISTANCE_REF_HPP

#include <cstdint>
#include <string>

std::uint32_t edit_distance_ref(
    const std::string& sequence_a,
    const std::string& sequence_b);

#endif