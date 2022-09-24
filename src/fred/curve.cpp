/*
Copyright 2020-2021 Dennis Rohde

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include <typeinfo>

#include "curve.hpp"
#include "simplification.hpp"

namespace fred {

Curve::Curve(const Points &points, const std::string &name) : Points(points), vstart{0}, vend{points.size() - 1}, name{name} {
    if (points.empty()) { 
        std::cerr << "warning: constructed empty curve" << std::endl;
        return; 
    }
    #if DEBUG
    std::cout << "constructed curve of complexity " << points.size() << std::endl;
    #endif
}


Curves Curves::simplify(const curve_size_t l, const bool approx = false) {
    Curves result(size(), l, Curves::dimensions());
    for (curve_number_t i = 0; i < size(); ++i) {
        if (approx) {
            Curve simplified_curve = Simplification::approximate_minimum_error_simplification(std::vector<Curve>::operator[](i), l);
            simplified_curve.set_name("Simplification of " + std::vector<Curve>::operator[](i).get_name());
            result[i] = simplified_curve;
        } else {
            Simplification::Subcurve_Shortcut_Graph graph(std::vector<Curve>::operator[](i));
            Curve simplified_curve = graph.minimum_error_simplification(l);
            simplified_curve.set_name("Simplification of " + std::vector<Curve>::operator[](i).get_name());
            result[i] = simplified_curve;
        }
        #if DEBUG
        std::cout << "Simplified curve " << i + 1 << "/" << size() << "." << std::endl;
        #endif
    }
    return result;
}

std::string Curve::repr() const {
    std::stringstream ss;
    ss << "fred.Curve '" << name << "' of complexity " << complexity() << " and " << dimensions() << " dimensions";
    return ss.str();
}

std::string Curves::repr() const {
    std::stringstream ss;
    ss << "fred.Curves collection with " << number() << " curves";
    return ss.str();
}

std::string Curve::str() const {
    std::stringstream ss;
    ss << name << std::endl;
    ss << *this;
    return ss.str();
}

std::string Curves::str() const {
    std::stringstream ss;
    ss << *this;
    return ss.str();
}

std::string Curve::get_name() const {
    return name;
}

void Curve::set_name(const std::string &name) {
    this->name = name;
}

std::ostream& operator<<(std::ostream &out, const Curve &curve) {
    if (curve.empty()) return out;
    out << "[";
    
    for (curve_size_t i = 0; i < curve.complexity() - 1; ++i) {
        out << curve[i] << ", ";
    }
    
    out << curve[curve.complexity() -1] << "]";

    return out;
}

std::ostream& operator<<(std::ostream &out, const Curves &curves) {
    if (curves.empty()) return out;
    out << "{";
    
    for (curve_number_t i = 0; i < curves.number() - 1; ++i) {
        out << curves[i] << ", ";
    }
    
    out << curves[curves.size() -1] << "}";

    return out;
}

}