
# How common-place are mediated effects?
yearly_mediate.data %>% View()
    ggplot(aes(x = year)) +
    geom_line(aes(y = mediate_count / general_count, colour = "blue"))
