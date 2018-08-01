
if(typeof(menu) == "string")
{
	document.write('\
			<ul>\
	')

	document.write('<li><a href="' + refdir + 'greens_functions.html" ')
	if(menu == 'greens_functions')
	{
		document.write('style="color: #b50404;')
	}
	document.write('">Greens functions</a> </li>')



	document.write('\
			</ul>\
	')
}

