
/*******************************************************************************
 * CGoGN: Combinatorial and Geometric modeling with Generic N-dimensional Maps  *
 * Copyright (C), IGG Group, ICube, University of Strasbourg, France            *
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

#ifndef CGOGN_RENDERING_TRANSFORM_FEEDBACKH_
#define CGOGN_RENDERING_TRANSFORM_FEEDBACK_H_

#include<vector>
#include<string>
#include<memory>

#include <cgogn/rendering/shaders/shader_program.h>



namespace cgogn
{

namespace rendering
{
template <typename SHADER>
class CGOGN_RENDERING_EXPORT TransformFeedback
{
	GLuint id_;
	typename SHADER::Param prg_param_;

	void internal_start(GLenum prim, const std::vector<VBO*>& vbos)
	{
		glEnable(GL_RASTERIZER_DISCARD);
		glBindTransformFeedback(GL_TRANSFORM_FEEDBACK, id_);
		for (GLuint i = 0; i < vbos.size(); ++i)
			glBindBufferBase(GL_TRANSFORM_FEEDBACK_BUFFER, i, vbos[i]->id());
		glBeginTransformFeedback(prim);
	}

public:
	TransformFeedback()
	{
		glCreateTransformFeedbacks(1,&id_);
		prg_param_ = SHADER::generate_param();
	}


	TransformFeedback(const TransformFeedback&) = delete;
	TransformFeedback& operator=(const TransformFeedback&) = delete;
	inline ~TransformFeedback() {}

	void start(GLenum prim, const std::vector<VBO*>& vbos)
	{
		prg_param_.bind();
		internal_start(prim,vbos);
	}

	void start(GLenum prim, const std::vector<VBO*>& vbos, const GLMat4& proj, const GLMat4& mv)
	{
		prg_param_.bind(proj,mv);
		internal_start(prim,vbos);
	}

	void stop()
	{
		prg_param_.release();
		glEndTransformFeedback();
		glBindTransformFeedback(GL_TRANSFORM_FEEDBACK,0);
		glDisable(GL_RASTERIZER_DISCARD);
	}
};

}
}
#endif
