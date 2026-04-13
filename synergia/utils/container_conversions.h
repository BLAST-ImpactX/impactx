#ifndef CONTAINER_CONVERSIONS_H_
#define CONTAINER_CONVERSIONS_H_

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace pybind11 {
    namespace detail {
        template <typename Type, typename Value>
        struct py_tuple_caster {
            using value_conv = make_caster<Value>;

            bool
            load(handle src, bool convert)
            {
                return false;
            }

            template <typename T>
            static handle
            cast(T&& src, return_value_policy policy, handle parent)
            {
                if (!std::is_lvalue_reference<T>::value)
                    policy =
                        return_value_policy_override<Value>::policy(policy);

                tuple t(src.size());
                size_t index = 0;

                for (auto&& value : src) {
                    auto value_ = reinterpret_steal<object>(value_conv::cast(
                        forward_like<T>(value), policy, parent));

                    if (!value_) return handle();

                    // steals a reference
                    PyTuple_SET_ITEM(
                        t.ptr(), (ssize_t)index++, value_.release().ptr());
                }

                return t.release();
            }

            PYBIND11_TYPE_CASTER(Type, _("Tuple[") + value_conv::name + _("]"));
        };
    }
}

#endif /* CONTAINER_CONVERSIONS_H_ */
