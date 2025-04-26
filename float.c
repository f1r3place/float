#include "return_codes.h"

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>

#define extract_mant(a, msize) ((a) & ((1 << msize) - 1))
#define extract_exp(a, msize, esize) (((a) >> (msize)) & ((1 << (esize)) - 1))
#define extract_sign(a, msize, esize) ((a) >> ((msize) + (esize)))

#define get_subnormal_mant(mant, exp) ((mant) << (-(exp) + 1))

enum cases
{
	normal,
	subnormal,
	nan,
	inf,
	zero
};

enum round
{
	towards_zero,
	ties_to_even,
	upward,
	downward
};

typedef struct
{
	uint8_t exp_size;
	uint8_t mant_size;
	uint8_t bias;
} precision_spec;

const precision_spec f16_spec = { .exp_size = 5, .mant_size = 10, .bias = 15 };
const precision_spec f32_spec = { .exp_size = 8, .mant_size = 23, .bias = 127 };

enum cases get_float_case(uint32_t input, precision_spec spec)
{
	uint8_t input_exp = extract_exp(input, spec.mant_size, spec.exp_size);
	uint32_t input_mant = extract_mant(input, spec.mant_size);
	if (input_exp == 0 && input_mant == 0)
	{
		return zero;
	}
	else if (input_exp == 0 && input_mant != 0)
	{
		return subnormal;
	}
	else if (input_exp == ((1 << spec.exp_size) - 1))
	{
		return input_mant == 0 ? inf : nan;
	}
	return normal;
}

uint32_t make_nan(precision_spec spec)
{
	uint32_t ret = 0;
	ret |= ((1 << spec.mant_size) - 2);
	ret |= ((1 << spec.exp_size) - 1) << spec.mant_size;
	return ret;
}

uint32_t make_inf(precision_spec spec, bool sign)
{
	uint32_t ret = 0;
	ret |= ((1 << spec.exp_size) - 1) << spec.mant_size;
	ret |= (sign << (spec.exp_size + spec.mant_size));
	return ret;
}

uint32_t make_max(precision_spec spec, bool sign)
{
	uint32_t ret = sign << (spec.mant_size + spec.exp_size);
	ret |= ((1 << spec.exp_size) - 2) << spec.mant_size;
	ret |= (1 << spec.mant_size) - 1;
	return ret;
}

uint32_t pack_float(precision_spec spec, bool sign, uint8_t exp, uint32_t mant)
{
	uint32_t ret = sign << (spec.mant_size + spec.exp_size) | exp << spec.mant_size | (mant & ((1 << spec.mant_size) - 1));
	return ret;
}

void float_invert_sign(uint32_t *ptr, precision_spec spec)
{
	*ptr ^= (1 << (spec.mant_size + spec.exp_size));
}

int16_t get_float_subnormal_exp(uint32_t input, precision_spec spec)
{
	int16_t initial_exp = -(spec.bias - 1);
	uint32_t input_mant = extract_mant(input, spec.mant_size);
	uint32_t initial_mant = input_mant;
	while ((initial_mant & (1 << spec.mant_size)) == 0)
	{
		initial_exp--;
		initial_mant <<= 1;
	}
	return initial_exp;
}

void subnormal_reroll(uint32_t a, int16_t *exp, uint64_t *mant, precision_spec spec)
{
	if (get_float_case(a, spec) == subnormal)
	{
		*exp = get_float_subnormal_exp(a, spec) + spec.bias;
		*mant = *mant << (-*exp + 1);
	}
}

void print_ufloat(uint32_t input, precision_spec spec)
{
	bool cur_sign = extract_sign(input, spec.mant_size, spec.exp_size);
	uint32_t cur_mant = extract_mant(input, spec.mant_size);
	uint8_t u_exp = extract_exp(input, spec.mant_size, spec.exp_size);
	cur_mant = cur_mant << (4 - (spec.mant_size % 4));
	int16_t cur_exp = u_exp - spec.bias;
	switch (get_float_case(input, spec))
	{
	case zero:
		if (cur_sign)
			printf("-");
		printf("0x0.");
		for (int i = 0; i < (spec.mant_size / 4) + 1; i++)
			printf("0");
		printf("p+0\n");
		break;
	case inf:
		if (cur_sign)
			printf("-");
		printf("inf\n");
		break;
	case nan:
		printf("nan\n");
		break;
	case subnormal:
		cur_exp = get_float_subnormal_exp(input, spec);
		cur_mant = get_subnormal_mant(cur_mant, cur_exp + spec.bias);
	case normal:
		if (cur_sign)
			printf("-");
		if (spec.mant_size == 10)
			printf("0x1.%x%x%x", (cur_mant >> 8) & 0xF, (cur_mant >> 4) & 0xF, cur_mant & 0xF);
		else
			printf(
				"0x1.%x%x%x%x%x%x",
				(cur_mant >> 20) & 0xF,
				(cur_mant >> 16) & 0xF,
				(cur_mant >> 12) & 0xF,
				(cur_mant >> 8) & 0xF,
				(cur_mant >> 4) & 0xF,
				cur_mant & 0xF);
		printf("p%+d\n", cur_exp);
		break;
	}
}

int error(int code, const char *message)
{
	fprintf(stderr, "%s\n", message);
	return code;
}

uint64_t float_round(uint64_t mant, bool sign, bool half_bit, bool sticky_bit, enum round round)
{
	bool last_bit = mant & 1;
	switch (round)
	{
	case towards_zero:
		break;
	case ties_to_even:
		mant += (sticky_bit && half_bit) || (half_bit && last_bit);
		break;
	case downward:
		mant += (sign && (sticky_bit | half_bit));
		break;
	case upward:
		mant += (!sign && (sticky_bit | half_bit));
		break;
	}
	return mant;
}

uint32_t float_exp_overflow(precision_spec spec, bool sign, enum round round)
{
	switch (round)
	{
	case towards_zero:
		return make_max(spec, sign);
	case upward:
		if (!sign)
		{
			return make_inf(spec, 0);
		}
		return make_max(spec, 1);
	case downward:
		if (!sign)
		{
			return make_max(spec, 0);
		}
		return make_inf(spec, 1);
	case ties_to_even:
		return make_inf(spec, sign);
	}
}

uint32_t float_multiply(uint32_t a, uint32_t b, enum round round, precision_spec spec)
{
	bool a_sign = extract_sign(a, spec.mant_size, spec.exp_size);
	bool b_sign = extract_sign(b, spec.mant_size, spec.exp_size);

	enum cases a_cases = get_float_case(a, spec);
	enum cases b_cases = get_float_case(b, spec);
	bool new_sign = a_sign ^ b_sign;
	if (a_cases == nan || b_cases == nan || ((a_cases == inf || b_cases == inf) && (b_cases == zero || a_cases == zero)))
	{
		return make_nan(spec);
	}
	else if (a_cases == inf || b_cases == inf)
	{
		return make_inf(spec, a_sign ^ b_sign);
	}
	else if (a_cases == zero || b_cases == zero)
	{
		return (a_sign ^ b_sign) << (spec.mant_size + spec.exp_size);
	}
	int16_t a_exp = extract_exp(a, spec.mant_size, spec.exp_size);
	int16_t b_exp = extract_exp(b, spec.mant_size, spec.exp_size);
	uint64_t a_mant = extract_mant(a, spec.mant_size), b_mant = extract_mant(b, spec.mant_size);
	if (a_cases == subnormal)
	{
		a_exp = get_float_subnormal_exp(a, spec) + spec.bias;
		a_mant = get_subnormal_mant(a_mant, a_exp);
	}

	if (b_cases == subnormal)
	{
		b_exp = get_float_subnormal_exp(b, spec) + spec.bias;
		b_mant = get_subnormal_mant(b_mant, b_exp);
	}

	int16_t new_exp = a_exp + b_exp - spec.bias;
	a_mant = a_mant | (1 << spec.mant_size);
	b_mant = b_mant | (1 << spec.mant_size);
	uint64_t unrounded_mant = a_mant * b_mant;
	uint8_t shift_amnt = 0;
	bool sticky_bit = 0;
	bool half_bit = 0;
	// this 'if' stays because the sticky_bit assignment condition is different
	// I know I could '&&' those two conditions but that would look REALLY messy
	if ((unrounded_mant & ((uint64_t)1 << (spec.mant_size * 2 + 1))))
	{
		new_exp++;
		sticky_bit = (unrounded_mant & ((1 << (spec.mant_size - 1)) - 1));
		unrounded_mant >>= 1;
	}
	sticky_bit |= ((unrounded_mant & ((1 << (spec.mant_size - 1)) - 1)) != 0);
	half_bit = (unrounded_mant & (1 << (spec.mant_size - 1))) != 0;
	unrounded_mant >>= spec.mant_size;
	uint64_t rounded_mant = (1 << spec.mant_size) ^ unrounded_mant;
	if (new_exp >= ((1 << spec.exp_size) - 1))
	{
		return float_exp_overflow(spec, new_sign, round);
	}
	if (new_exp <= 0)
	{
		shift_amnt -= new_exp - 1;
		if (shift_amnt <= spec.mant_size)
		{
			sticky_bit |= half_bit;
			half_bit = (unrounded_mant & (1 << (shift_amnt - 1))) != 0;
		}
		new_exp = 0;
		rounded_mant |= 1 << spec.mant_size;
		if (shift_amnt <= (spec.mant_size + 1))
		{
			sticky_bit = ((((1 << (shift_amnt - 1)) - 1) & unrounded_mant) != 0) | sticky_bit;
			rounded_mant >>= shift_amnt;
			half_bit |= (shift_amnt == spec.mant_size + 1);
		}
		else
		{
			sticky_bit = 1;
			half_bit = 0;
			rounded_mant = 0;
		}
	}
	rounded_mant = float_round(rounded_mant, new_sign, half_bit, sticky_bit, round);
	new_exp += (rounded_mant & (1 << spec.mant_size)) != 0;
	return pack_float(spec, new_sign, new_exp, rounded_mant);
}

uint32_t float_subtract(uint32_t a, uint32_t b, enum round round, precision_spec spec)
{
	bool a_sign = extract_sign(a, spec.mant_size, spec.exp_size);
	bool b_sign = extract_sign(b, spec.mant_size, spec.exp_size);

	enum cases a_cases = get_float_case(a, spec);
	enum cases b_cases = get_float_case(b, spec);

	if (a_cases == zero)
	{
		b_sign = (a_sign ^ b_sign);
		return pack_float(spec, b_sign, extract_exp(a, spec.mant_size, spec.exp_size), extract_mant(a, spec.mant_size));
	}
	else if (b_cases == zero)
	{
		return a;
	}

	if (a_cases == nan || b_cases == nan || (a_cases == inf && b_cases == inf && (a_sign == b_sign)))
	{
		return make_nan(spec);
	}

	if (a_cases == inf)
	{
		return make_inf(spec, a_sign);
	}

	if (b_cases == inf)
	{
		return make_inf(spec, !b_sign);
	}

	int16_t a_exp = extract_exp(a, spec.mant_size, spec.exp_size);
	int16_t b_exp = extract_exp(b, spec.mant_size, spec.exp_size);
	uint64_t a_mant = extract_mant(a, spec.mant_size), b_mant = extract_mant(b, spec.mant_size);
	if (a_cases == subnormal)
	{
		a_exp = get_float_subnormal_exp(a, spec) + spec.bias;
		a_mant = get_subnormal_mant(a_mant, a_exp);
	}

	if (b_cases == subnormal)
	{
		b_exp = get_float_subnormal_exp(b, spec) + spec.bias;
		b_mant = get_subnormal_mant(b_mant, b_exp);
	}

	bool new_sign = a_sign;
	if (a_exp < b_exp || (a_exp == b_exp && a_mant < b_mant))
	{
		uint32_t t = b;
		b = a;
		a = t;
		int16_t t_exp = b_exp;
		b_exp = a_exp;
		a_exp = t_exp;
		uint64_t t_mant = b_mant;
		b_mant = a_mant;
		a_mant = t_mant;
		new_sign = !new_sign;
	}	 // after this we can safely assume a >= b

	uint16_t exp_diff = a_exp - b_exp;
	int16_t new_exp = a_exp;
	a_mant = a_mant | 1 << (spec.mant_size);
	b_mant = b_mant | 1 << (spec.mant_size);
	uint64_t unrounded_mant;
	bool sticky_bit = 0, half_bit = 0;
	if (exp_diff <= (spec.mant_size + 1))
	{
		a_mant <<= exp_diff;
		unrounded_mant = a_mant - b_mant;
		if (unrounded_mant == 0)
		{
			if (round == downward)
				return 1 << (spec.exp_size + spec.mant_size);
			return 0;
		}
		while (!(unrounded_mant & ((uint64_t)1 << (spec.mant_size + exp_diff))))
		{
			new_exp--;
			unrounded_mant <<= 1;
		}
		if (exp_diff > 0)
		{
			half_bit = ((1 << (exp_diff - 1)) & unrounded_mant) != 0;
		}
		if (exp_diff > 1)
		{
			sticky_bit = (((1 << (exp_diff - 1)) - 1) & unrounded_mant) != 0;
		}
		unrounded_mant >>= exp_diff;
	}
	else
	{
		unrounded_mant = a_mant;
		unrounded_mant -= 1;
		sticky_bit = 1;
		half_bit = 1;
		if (!(unrounded_mant & (1 << spec.mant_size)))
		{
			new_exp--;
			unrounded_mant <<= 1;
			unrounded_mant |= 1;
			half_bit &= !(exp_diff == (spec.mant_size + 2) && b_mant != (1 << spec.mant_size));
		}
	}
	uint64_t rounded_mant = float_round(unrounded_mant, new_sign, half_bit, sticky_bit, round);
	new_exp += (!(rounded_mant & (1 << spec.mant_size)));
	unrounded_mant >>= (!(rounded_mant & (1 << spec.mant_size)));
	if (new_exp <= 0)
	{
		if ((-new_exp + 1) <= spec.mant_size)
			rounded_mant >>= (-new_exp + 1);
		else
			rounded_mant = 0;
		new_exp = 0;
	}
	return pack_float(spec, new_sign, new_exp, rounded_mant);
}

uint32_t float_addition(uint32_t a, uint32_t b, enum round round, precision_spec spec)
{
	bool a_sign = extract_sign(a, spec.mant_size, spec.exp_size);
	bool b_sign = extract_sign(b, spec.mant_size, spec.exp_size);

	enum cases a_cases = get_float_case(a, spec);
	enum cases b_cases = get_float_case(b, spec);

	if (a_cases == zero)
	{
		return b;
	}
	else if (b_cases == zero)
	{
		return a;
	}

	if (a_cases == nan || b_cases == nan || (a_cases == inf && b_cases == inf && (a_sign != b_sign)))
	{
		return make_nan(spec);
	}

	if (a_cases == inf)
	{
		return make_inf(spec, a_sign);
	}

	if (b_cases == inf)
	{
		return make_inf(spec, b_sign);
	}

	int16_t a_exp = extract_exp(a, spec.mant_size, spec.exp_size);
	int16_t b_exp = extract_exp(b, spec.mant_size, spec.exp_size);
	uint64_t a_mant = extract_mant(a, spec.mant_size), b_mant = extract_mant(b, spec.mant_size);
	if (a_cases == subnormal)
	{
		a_exp = get_float_subnormal_exp(a, spec) + spec.bias;
		a_mant = get_subnormal_mant(a_mant, a_exp);
	}

	if (b_cases == subnormal)
	{
		b_exp = get_float_subnormal_exp(b, spec) + spec.bias;
		b_mant = get_subnormal_mant(b_mant, b_exp);
	}

	if (a_exp < b_exp)
	{
		uint32_t t = b;
		b = a;
		a = t;
		int16_t t_exp = b_exp;
		b_exp = a_exp;
		a_exp = t_exp;
		uint64_t t_mant = b_mant;
		b_mant = a_mant;
		a_mant = t_mant;
	}	 // after this we can safely assume a.exp >= b.exp

	uint16_t exp_diff = a_exp - b_exp;	  // >= 0
	int16_t new_exp = a_exp;
	a_mant = a_mant | (1 << spec.mant_size);
	b_mant = b_mant | (1 << spec.mant_size);

	bool sticky_bit = 0;
	bool half_bit = 0;
	if (exp_diff > 0)
	{
		half_bit = ((1 << (exp_diff - 1)) & b_mant) != 0;
		sticky_bit = (((1 << (exp_diff - 1)) - 1) & b_mant) != 0 && (exp_diff > 1);
	}
	if (exp_diff == (spec.mant_size + 1))
	{
		half_bit = 1;
		sticky_bit = (((1 << spec.mant_size) - 1) & b_mant) != 0;
	}
	if (exp_diff > spec.mant_size + 1)
	{
		half_bit = 0;
		sticky_bit = 1;
	}
	if (exp_diff > spec.mant_size)
		b_mant = 0;
	else
		b_mant >>= exp_diff;
	uint64_t unrounded_mant = a_mant + b_mant;
	// this 'if' stays because of sticky_bit and half_bit assignment
	// possible replacement without 'if' looks like:
	//	 sticky_bit |= half_bit && ((unrounded_mant & (1 << (spec.mant_size + 1))) != 0);
	// but it's very hard to read
	// even if the condition will be in a separate variable
	//	 bool x = ((unrounded_mant & (1 << ...));
	//	 new_exp += x;
	//	 sticky_bit |= x && half_bit;
	//	 half_bit = (x && (unrounded_mant & 1)) || (!x && half_bit);
	//	 unrounded_mant >> x;
	// this is only two lines less but is a lot less readable
	if ((unrounded_mant & (1 << (spec.mant_size + 1))))
	{
		new_exp++;
		sticky_bit = sticky_bit | half_bit;
		half_bit = unrounded_mant & 1;
		unrounded_mant >>= 1;
	}
	uint64_t rounded_mant = float_round(unrounded_mant, a_sign, half_bit, sticky_bit, round);
	new_exp += (rounded_mant & (1 << (spec.mant_size + 1))) != 0;
	rounded_mant >>= (rounded_mant & (1 << (spec.mant_size + 1))) != 0;
	if (new_exp <= 0)
	{
		if (-new_exp + 1 <= spec.mant_size)
			rounded_mant >>= (-new_exp + 1);
		else
			rounded_mant = 0;
		new_exp = 0;
	}
	if (new_exp >= ((1 << spec.exp_size) - 1))
	{
		return float_exp_overflow(spec, a_sign, round);
	}
	return pack_float(spec, a_sign, new_exp, rounded_mant);
}

uint32_t float_add_diff(uint32_t a, uint32_t b, enum round round, precision_spec spec)
{
	bool a_sign = extract_sign(a, spec.mant_size, spec.exp_size);
	bool b_sign = extract_sign(b, spec.mant_size, spec.exp_size);
	if (a_sign == b_sign)
	{
		return float_addition(a, b, round, spec);
	}
	if (get_float_case(a, spec) == zero && get_float_case(b, spec) == zero)
	{
		return (1 << (spec.mant_size + spec.exp_size)) * (round == downward);
	}
	if (get_float_case(a, spec) == zero)
	{
		return b;
	}
	else if (get_float_case(b, spec) == zero)
	{
		return a;
	}
	if (b_sign)
	{
		float_invert_sign(&b, spec);
		return float_subtract(a, b, round, spec);
	}
	else
	{
		float_invert_sign(&a, spec);
		return float_subtract(b, a, round, spec);
	}
}

uint32_t float_sub_diff(uint32_t a, uint32_t b, enum round round, precision_spec spec)
{
	float_invert_sign(&b, spec);
	return float_add_diff(a, b, round, spec);
}

uint32_t float_divide(uint32_t a, uint32_t b, enum round round, precision_spec spec)
{
	bool a_sign = extract_sign(a, spec.mant_size, spec.exp_size);
	bool b_sign = extract_sign(b, spec.mant_size, spec.exp_size);

	enum cases a_cases = get_float_case(a, spec);
	enum cases b_cases = get_float_case(b, spec);

	if (a_cases == nan || b_cases == nan || (a_cases == zero && b_cases == zero) || (a_cases == inf && b_cases == inf))
	{
		return make_nan(spec);
	}

	if (b_cases == zero || a_cases == inf)
	{
		return make_inf(spec, a_sign ^ b_sign);
	}

	if (a_cases == zero || b_cases == inf)
	{
		return ((a_sign ^ b_sign) << (spec.exp_size + spec.mant_size));
	}

	int16_t a_exp = extract_exp(a, spec.mant_size, spec.exp_size);
	int16_t b_exp = extract_exp(b, spec.mant_size, spec.exp_size);
	uint64_t a_mant = extract_mant(a, spec.mant_size), b_mant = extract_mant(b, spec.mant_size);
	if (a_cases == subnormal)
	{
		a_exp = get_float_subnormal_exp(a, spec) + spec.bias;
		a_mant = get_subnormal_mant(a_mant, a_exp);
	}

	if (b_cases == subnormal)
	{
		b_exp = get_float_subnormal_exp(b, spec) + spec.bias;
		b_mant = get_subnormal_mant(b_mant, b_exp);
	}

	bool new_sign = a_sign ^ b_sign;
	int16_t new_exp = a_exp - b_exp + spec.bias;
	a_mant = a_mant | 1 << spec.mant_size;
	b_mant = b_mant | 1 << spec.mant_size;
	// we need exactly two bits of additional precision (to know half_bit after
	// normalization)
	a_mant <<= (spec.mant_size + 3);
	bool sticky_bit, half_bit;
	uint64_t unrounded_mant = a_mant / b_mant;
	new_exp -= ((unrounded_mant & (1 << (spec.mant_size + 3))) == 0);
	unrounded_mant <<= ((unrounded_mant & (1 << (spec.mant_size + 3))) == 0);
	half_bit = (unrounded_mant & (1 << 2)) != 0;
	sticky_bit = (a_mant % b_mant) != 0;
	unrounded_mant >>= 3;
	if (new_exp >= ((1 << spec.exp_size) - 1))
	{
		return float_exp_overflow(spec, new_sign, round);
	}
	if (new_exp <= 0)
	{
		uint16_t shift_amnt = (-new_exp + 1);
		if (shift_amnt <= (spec.mant_size + 1))
		{
			sticky_bit = half_bit | sticky_bit | (unrounded_mant & ((1 << (shift_amnt - 1)) - 1));
			half_bit = (unrounded_mant & (1 << (shift_amnt - 1))) != 0 || (shift_amnt == spec.mant_size + 1);
		}
		else
		{
			sticky_bit = 1;
			half_bit = 0;
		}
		if (shift_amnt <= spec.mant_size)
			unrounded_mant >>= shift_amnt;
		else
			unrounded_mant = 0;
		new_exp = 0;
		uint32_t norm_mant = float_round(unrounded_mant, new_sign, half_bit, sticky_bit, round);
		new_exp += ((norm_mant & (1 << spec.mant_size)) != 0);
	}
	uint32_t rounded_mant = float_round(unrounded_mant, new_sign, half_bit, sticky_bit, round);
	return pack_float(spec, new_sign, new_exp, rounded_mant);
}

int user_interface(int argv, char **argc)
{
	if (argv == 1)
	{
		printf("usage: float <precision (h | f)> <rounding> <first number> [<operation> <second number>]\n");
		return SUCCESS;
	}
	else if (argv == 4 || argv == 6)
	{
		enum round round;
		sscanf(argc[2], "%u", &round);
		if (round < 0 || round > 3)
		{
			return error(ERROR_ARGUMENTS_INVALID, "invalid rounding");
		}
		char c;
		c = *argc[1];
		if (*(argc[1] + 1) != '\0')
			return error(ERROR_ARGUMENTS_INVALID, "invalid precision type");
		if (c != 'h' && c != 'f')
		{
			return error(ERROR_ARGUMENTS_INVALID, "invalid precision type");
		}
		uint32_t n1;
		sscanf(argc[3], "0x%x", &n1);
		precision_spec passed_spec = c == 'h' ? f16_spec : f32_spec;
		if (argv == 4)
		{
			print_ufloat(n1, passed_spec);
			return SUCCESS;
		}

		char op;
		op = *argc[4];
		if (*(argc[4] + 1) != '\0')
		{
			return error(ERROR_ARGUMENTS_INVALID, "invalid operation");
		}
		uint32_t n2;
		sscanf(argc[5], "0x%x", &n2);
		uint32_t res;
		switch (op)
		{
		case '*':
			res = float_multiply(n1, n2, round, passed_spec);
			break;
		case '+':
			res = float_add_diff(n1, n2, round, passed_spec);
			break;
		case '-':
			res = float_sub_diff(n1, n2, round, passed_spec);
			break;
		case '/':
			res = float_divide(n1, n2, round, passed_spec);
			break;
		default:
			return error(ERROR_ARGUMENTS_INVALID, "invalid operation");
		}
		print_ufloat(res, passed_spec);
	}
	else
	{
		return error(ERROR_ARGUMENTS_INVALID, "incorrect argument amount");
	}
	return SUCCESS;
}

int main(int argv, char **argc)
{
	return user_interface(argv, argc);
}
