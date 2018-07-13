
if(typeof(menu) == "string")
{
	document.write('\
			<ul>\
	')


	document.write('<li><a href="' + refdir + 'mapping.html" ')
	if(menu == 'mapping')
	{
		document.write('style="color: #b50404;')
	}
	document.write('">Mapping</a> </li>')


	document.write('\
			</ul>\
	')
}


