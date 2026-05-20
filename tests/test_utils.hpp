#ifndef JABBER_TEST_UTILS
#define JABBER_TEST_UTILS

#include <catch2/generators/catch_generators_all.hpp>
namespace jabber_test
{

/// Concept for enums of type std::uint8_t with a ::Size.
template<typename T>
concept OptionEnum =
   std::is_enum_v<T> &&
   std::same_as<std::underlying_type_t<T>, std::uint8_t> &&
   requires { T::Size; };

/// Custom generator for \ref OptionEnum.
template<OptionEnum T>
class OptionGenerator : public Catch::Generators::IGenerator<T>
{
private:
   T enumerator_val_;
   std::uint8_t enum_idx_ = 0;
public:
   bool next() override
   {
      if (++enum_idx_ < static_cast<std::uint8_t>(T::Size))
      {
         enumerator_val_ = static_cast<T>(enum_idx_);
         return true;
      }
      else
      {
         return false;
      }
   };

   T const &get() const override
   {
      return enumerator_val_;
   }
};

/**
 * @brief Generator-wrapper creator for list over all options by
 *  \ref OptionEnum. Use by \c GENERATE(options<T>()).
 *
 */
template <OptionEnum T>
Catch::Generators::GeneratorWrapper<T> options()
{
   return Catch::Generators::GeneratorWrapper<T>(
             Catch::Detail::make_unique<OptionGenerator<T>>());
}

/// Random generator for a \ref OptionEnum enumerator.
template<OptionEnum T>
class RandomOptionGenerator : public Catch::Generators::IGenerator<T>
{
private:
   Catch::Generators::GeneratorWrapper<std::uint8_t> ri_gen_;
   T option_;
public:
   RandomOptionGenerator()
      : ri_gen_(Catch::Generators::random<std::uint8_t>(
                   0,static_cast<std::uint8_t>(T::Size)-1))
   {
      static_cast<void>(next());
   }

   T const &get() const override
   {
      return option_;
   }

   bool next() override
   {
      ri_gen_.next();
      option_ = static_cast<T>(ri_gen_.get());
      return true;
   }
};

/**
 * @brief Generator-wrapper creator for random option of
 * \ref OptionEnum. Use by \c random_option<E>().
 */
template<OptionEnum T>
Catch::Generators::GeneratorWrapper<T> random_option()
{
   return Catch::Generators::GeneratorWrapper<T>(
             Catch::Detail::make_unique<RandomOptionGenerator<T>>());
}

} // namespace jabber_test

#endif // JABBER_TEST_UTILS
