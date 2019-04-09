/*******************************************************************************
* CGoGN: Combinatorial and Geometric modeling with Generic N-dimensional Maps  *
* Copyright (C) 2015, IGG Group, ICube, University of Strasbourg, France       *
*                                                                              *
* This library is free software; you can redistribute it and/or modify it      *
* under the terms of the GNU Lesser General Public License as published by the *
* Free Software Foundation; either version 2.1 of the License, or (at your     *
* option) any later version.                                                   *
*                                                                              *
* This library is distributed in the hope that it will be useful, but WITHOUT  *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or        *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License  *
* for more details.                                                            *
*                                                                              *
* You should have received a copy of the GNU Lesser General Public License     *
* along with this library; if not, write to the Free Software Foundation,      *
* Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.           *
*                                                                              *
* Web site: http://cgogn.unistra.fr/                                           *
* Contact information: cgogn@unistra.fr                                        *
*                                                                              *
*******************************************************************************/

#ifndef CGOGN_CORE_UTILS_ASSERT_H_
#define CGOGN_CORE_UTILS_ASSERT_H_

#include <cgogn/core/cgogn_core_export.h>

#include <cgogn/core/utils/definitions.h>

#include <string>

#if defined (WIN32) && !defined(__func__)
#define __func__ __FUNCTION__
#endif

/**
 * \file cgogn/core/basic/assert.h
 * \brief Assertion checking mechanism.
 *
 * Allows the user to add a specific message to output
 */
namespace cgogn
{

/**
 * Prints an assertion failure.
 * This function is called when a boolean condition is not met.
 * It prints an error message and terminates the program.
 * \param[in] expression string representation of the condition.
 * \param[in] message string information message to print out.
 * \param[in] file_name file where the assertion failed.
 * \param[in] function_name function where the assertion failed.
 * \param[in] line_number line where the assertion failed.
 */
CGOGN_NORETURN CGOGN_CORE_EXPORT void assertion_failed(
	const std::string& expression,
	const std::string& message,
	const std::string& file_name,
	const std::string& function_name,
	int line_number
);

/**
 * Prints an unreachable location failure.
 * This function is called when execution reaches a point that
 * should not be reached. It prints an error message and
 * terminates the program.
 * \param[in] message string information message to print out.
 * \param[in] file_name file where the assertion failed.
 * \param[in] function_name function where the assertion failed.
 * \param[in] line_number line where the assertion failed.
 */
CGOGN_NORETURN CGOGN_CORE_EXPORT void should_not_have_reached(
	const std::string& message,
	const std::string& file_name,
	const std::string& function_name,
	int line_number
);

} // namespace cgogn

/**
 * \brief Verifies that a condition is met.
 * \details
 * \param[in] x the boolean expression of the condition
 * \see assertion_failed()
 */
#define _internal_cgogn_assert(x) 										\
do { 																	\
	if(!(x)) 															\
	{			 														\
		cgogn::assertion_failed(#x, "", __FILE__, __func__, __LINE__);	\
	} 																	\
} while (false)

/**
 * \brief Verifies that a condition is met and take a specific message.
 * \details
 * \param[in] x the boolean expression of the condition
 * \param[in] msg the specific message about the condition
 * \see assertion_failed()
 */
#define _internal_cgogn_message_assert(x, msg)							\
do { 																	\
	if(!(x)) 															\
	{			 														\
		cgogn::assertion_failed(#x, msg, __FILE__, __func__, __LINE__);	\
	} 																	\
} while (false)

/**
 * \brief Sets a non reachable point in the program
 * \details
 * \param[in] msg the specific information message
 */
#define _internal_cgogn_assert_not_reached(msg)							\
do {																	\
	cgogn::should_not_have_reached(msg, __FILE__, __func__, __LINE__);	\
} while (false)

/**
 * \brief Verifies that the required contract condition is met.
 * \details
 *
 * \param[in] x the boolean expression of the condition
 * \see assertion_failed()
 */
#define _internal_cgogn_require(x) 										\
do { 																	\
	if(!(x)) 															\
	{			 														\
		cgogn::assertion_failed(#x, "", __FILE__, __func__, __LINE__);	\
	} 																	\
} while (false)

/**
 * \brief Verifies that the ensured contract condition is met.
 * \details
 *
 * \param[in] x the boolean expression of the condition
 * \see assertion_failed()
 */
#define _internal_cgogn_ensure(x) 										\
do { 																	\
	if(!(x)) 															\
	{			 														\
		cgogn::assertion_failed(#x, "", __FILE__, __func__, __LINE__);	\
	} 																	\
} while (false)

/**
 * \brief Verifies that the invariant contract condition is met.
 * \details
 *
 * \param[in] x the boolean expression of the condition
 * \see assertion_failed()
 */
#define _internal_cgogn_invariant(x) 									\
do { 																	\
	if(!(x)) 															\
	{			 														\
		cgogn::assertion_failed(#x, "", __FILE__, __func__, __LINE__);	\
	} 																	\
} while (false)

/**
 * \def cgogn_assert(x)
 * \copydoc cgogn_assert()
 * \note This assertion check is only active in debug mode.
 */
 /**
 * \def cgogn_message_assert(x, msg)
 * \copydoc cgogn_assert()
 * \note This assertion check is only active in debug mode.
 */
#ifdef CGOGN_DEBUG
	#define cgogn_assert(x) _internal_cgogn_assert(x)
	#define cgogn_message_assert(x, msg) _internal_cgogn_message_assert(x, msg)
	#define cgogn_assert_not_reached(msg) _internal_cgogn_assert_not_reached(msg)
	#define cgogn_assert_not_reached_false(msg) _internal_cgogn_assert_not_reached(msg); return false;
	#define cgogn_require(x) _internal_cgogn_require(x)
	#define cgogn_ensure(x) _internal_cgogn_ensure(x)
	#define cgogn_invariant(x) _internal_cgogn_invariant(x)
#else
	#define cgogn_assert(x)
	#define cgogn_message_assert(x, msg)
	#define cgogn_assert_not_reached(msg)
	#define cgogn_assert_not_reached_false(msg) return false;
	#define cgogn_require(x)
	#define cgogn_ensure(x)
	#define cgogn_invariant(x)
#endif

/**
 * \def parano_assert(x)
 * \copydoc cgogn_assert()
 * \note This assertion check is only active in parano mode.
 */
 /**
 * \def parano_message_assert(x, msg)
 * \copydoc cgogn_assert()
 * \note This assertion check is only active in parano mode.
 */
#ifdef CGOGN_PARANO
	#define parano_assert(x) _internal_cgogn_assert(x)
	#define parano_message_assert(x, msg) _internal_cgogn_message_assert(x, msg)
#else
	#define parano_assert(x)
	#define parano_message_assert(x, msg)
#endif

#endif // CGOGN_CORE_UTILS_ASSERT_H_
