
if(typeof(menu) == "string")
{
	document.write('\
			<ul>\
	')

	document.write('<li><a href="' + refdir + 'map2d.html" ')
	if(menu == 'map2d')
	{
		document.write('style="color: #b50404;')
	}
	document.write('">map2d</a> </li>')


	document.write('\
			</ul>\
	')
}


