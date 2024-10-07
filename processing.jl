using ImageMagick
using FileIO
using Cairo, Rsvg

function svg_to_png(input_file::AbstractString, output_file::AbstractString)
    r = Rsvg.handle_new_from_file(input_file)
    d = Rsvg.handle_get_dimensions(r)
    cs = Cairo.CairoImageSurface(d.width, d.height, Cairo.FORMAT_ARGB32)
    c = Cairo.CairoContext(cs)
    Rsvg.handle_render_cairo(c, r)
    Cairo.write_to_png(cs, output_file)

end

svg_to_png("imgs/Course_Icon_LinAlgebra.svg", "imgs/math557_logo.png")