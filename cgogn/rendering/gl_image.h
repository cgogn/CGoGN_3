
#ifndef CGOGN_RENDERING_GL_IMAGE_H_
#define CGOGN_RENDERING_GL_IMAGE_H_

#include <cgogn/rendering/cgogn_rendering_export.h>
#include <cgogn/rendering/types.h>
#include <vector>

namespace cgogn
{

namespace rendering
{

class CGOGN_RENDERING_EXPORT GLImage
{
	uint8* data_;
	int32 width_;
	int32 height_;
	int32 channels_; // 1: grey, 3: RGB, 4: RGBA
	bool stb_loaded_;

public:
	GLImage(int32 width, int32 height, int32 channels);
	GLImage(const std::string& filename);
	~GLImage();

	CGOGN_NOT_COPYABLE_NOR_MOVABLE(GLImage);

	void save(const std::string& filename, bool flip_y = false) const;

	inline int32 width() const
	{
		return width_;
	}
	inline int32 height() const
	{
		return height_;
	}
	inline int32 depth() const
	{
		return channels_;
	}
	inline const uint8* data() const
	{
		return data_;
	}
	
	inline void copy_pixels_data(uint8* ptr)
	{
		std::memcpy(data_, ptr, width_ * height_ * channels_);
	}

	inline void copy_pixels_data(const std::vector<uint8> pix)
	{
		auto sz = width_ * height_ * channels_;
		if (pix.size() >= sz)
			std::memcpy(data_, pix.data(), sz);
	}

	inline void set_pixel(int32 x, int32 y, const GLColor& col)
	{
		uint8* ptr = data_ + (uint64(channels_) * (y * width_ + x));
		for (int32 i = 0; i < channels_; ++i)
			*ptr++ = uint8(255 * col[i]);
	}
};

} // namespace rendering

} // namespace cgogn

#endif // CGOGN_RENDERING_GL_IMAGE_H_
